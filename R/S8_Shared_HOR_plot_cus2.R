#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  # ggplot2 provides scale_*_viridis_* via viridisLite
})

options(scipen = 999)
`%!in%` <- Negate(`%in%`)

# ---------- CLI ----------
args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "$|^", flag, "="), args)
  if (!length(hit)) return(default)
  a <- args[hit[1]]
  if (grepl("=", a)) sub(paste0("^", flag, "="), "", a) else {
    i <- hit[1] + 1; if (i <= length(args)) args[i] else default
  }
}

OUT_DIR <- argval("--outdir", default = NULL)
if (is.null(OUT_DIR) || OUT_DIR == "") stop("Missing --outdir")
OUT_DIR <- normalizePath(OUT_DIR, mustWork = TRUE)

PLOTV  <- toupper(argval("--plotv", default = ""))  # "", V1, V2, V3
if (nzchar(PLOTV) && PLOTV %!in% c("V1","V2","V3")) {
  message("WARN: --plotv should be one of V1,V2,V3. Got: ", PLOTV,
          ". Proceeding with default (V1).")
  PLOTV <- ""
}

STRING_DIR <- file.path(OUT_DIR, "string_share_HOR")
if (!dir.exists(STRING_DIR)) stop("Not found: ", STRING_DIR)

# ---------- read data (support both file name styles) ----------
data_files <- dir(STRING_DIR,
                  pattern = "_HOR_patterns_bed\\.csv$",
                  full.names = TRUE)
if (!length(data_files)) {
  message("No *_HOR_patterns_bed.csv found under ", STRING_DIR)
  quit(save = "no", status = 0)
}
data_files <- sort(data_files)

# attempt to extract threshold from common patterns:
#   .../chrN_93_bins_nam_groups_HOR_patterns_bed.csv
#   .../_93_bins_..._HOR_patterns_bed.csv (no chr in front)
extract_thr <- function(p) {
  b <- basename(p)
  # try chrN_XX_bins first
  m1 <- regmatches(b, regexec("chr\\w+_(\\d+)_bins", b))[[1]]
  if (length(m1) >= 2) return(as.integer(m1[2]))
  # fallback: _XX_bins anywhere
  m2 <- regmatches(b, regexec("_(\\d+)_bins", b))[[1]]
  if (length(m2) >= 2) return(as.integer(m2[2]))
  NA_integer_
}

normalize_types <- function(d) {
  d %>%
    mutate(
      chr           = as.character(chr),
      array_start   = suppressWarnings(as.numeric(array_start)),
      array_end     = suppressWarnings(as.numeric(array_end)),
      bin_start     = suppressWarnings(as.numeric(bin_start)),
      bin_end       = suppressWarnings(as.numeric(bin_end)),
      pattern_start = suppressWarnings(as.numeric(pattern_start)),
      pattern_end   = suppressWarnings(as.numeric(pattern_end))
    )
}

dat_list <- lapply(seq_along(data_files), function(i) {
  thr <- extract_thr(data_files[i])
  d <- tryCatch(read.csv(data_files[i], stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- normalize_types(d)
  if (!nrow(d)) return(NULL)
  d$threshold <- thr
  d
})
dat_list <- Filter(Negate(is.null), dat_list)
if (!length(dat_list)) {
  message("All *_HOR_patterns_bed.csv are empty. Nothing to plot.")
  quit(save = "no", status = 0)
}

df <- dplyr::bind_rows(dat_list)

# ---------- (optional) keep patterns spanning >1 bin per chr/threshold ----------
# unique by chr, bin_start, HOR_pattern, threshold; require n > 1
span_filt <- df %>%
  select(chr, bin_start, HOR_pattern, threshold) %>%
  distinct() %>%
  group_by(chr, HOR_pattern, threshold) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

df_merge <- df %>%
  inner_join(span_filt, by = c("chr","HOR_pattern","threshold")) %>%
  # clean up columns (keep original fields; drop the temp 'n')
  select(-n)

# if the filter removed everything, fall back to original df
if (!nrow(df_merge)) df_merge <- df

# ---------- derive fields ----------
df_all <- df_merge %>%
  mutate(threshold = as.integer(threshold)) %>%
  arrange(threshold, HOR_pattern) %>%
  group_by(threshold, HOR_pattern) %>%
  mutate(order_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(
    array_start   = as.numeric(array_start),
    array_end     = as.numeric(array_end),
    pattern_start = as.numeric(pattern_start),
    pattern_end   = as.numeric(pattern_end),
    pos           = (pattern_start + pattern_end) / 2
  )

if (!nrow(df_all)) {
  message("No rows to plot after combining inputs. Exiting.")
  quit(save = "no", status = 0)
}

df_all$HOR_pattern <- str_replace_all(df_all$HOR_pattern, "\\|", " ")
df_all <- df_all %>% mutate(mon_size = stringr::str_count(HOR_pattern, "\\S+"))

lab_map <- setNames(df_all$HOR_pattern, df_all$order_id)
xlim_all <- c(min(df_all$array_start, na.rm = TRUE),
              max(df_all$array_end,   na.rm = TRUE))

# pick a chromosome label (if single) for filenames
uniq_chr <- unique(na.omit(df_all$chr))
chr_tag  <- if (length(uniq_chr) == 1) uniq_chr else "all"

# thresholds palette
thresholds <- 90:99
thresh_cols <- c(
  "90"="#5AB4AC","91"="#f2a287","92"="#aa00ff","93"="#008e17","94"="#F0027F",
  "95"="#bc9000","96"="maroon","97"="#3c3c3c","98"="#0099CC","99"="red"
)

# ---------- PLOTS ----------
# V1: bin_start on x, manual threshold palette
p <- ggplot(df_all, aes(x = bin_start, y = order_id, color = factor(threshold))) +
  geom_point(alpha = 0.85, size = 1) +
  scale_y_continuous(
    breaks = as.numeric(names(lab_map)),
    labels = lab_map,
    limits = c(0, max(df_all$order_id, na.rm = TRUE) + 1),
    expand = expansion(mult = c(0, 0), add = 0)
  ) +
  scale_x_continuous(limits = xlim_all) +
  scale_color_manual(values = thresh_cols,
                     breaks = thresholds,
                     name = "Threshold") +
  labs(x = "Genomic position", y = "HOR Pattern") +
  theme_bw(base_size = 40) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 40),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.ticks.length.y = unit(1, "mm")
  )

# V2: bin_start on x, discrete viridis palette
p2 <- ggplot(df_all, aes(x = bin_start, y = order_id, color = factor(threshold))) +
  geom_point(alpha = 0.85, size = 1) +
  scale_y_continuous(
    breaks = as.numeric(names(lab_map)),
    labels = lab_map,
    limits = c(0, max(df_all$order_id, na.rm = TRUE) + 1),
    expand = expansion(mult = c(0, 0), add = 0)
  ) +
  scale_x_continuous(limits = xlim_all) +
  scale_color_viridis_d(name = "Threshold", option = "plasma") +
  labs(x = "Genomic position", y = "HOR Pattern") +
  theme_bw(base_size = 40) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 40),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.ticks.length.y = unit(1, "mm")
  )

# V3: pattern center on x, discrete viridis palette
p3 <- ggplot(df_all, aes(x = (pattern_start + pattern_end)/2, y = order_id, color = factor(threshold))) +
  geom_point(alpha = 0.85, size = 1) +
  scale_y_continuous(
    breaks = as.numeric(names(lab_map)),
    labels = lab_map,
    limits = c(0, max(df_all$order_id, na.rm = TRUE) + 1),
    expand = expansion(mult = c(0, 0), add = 0)
  ) +
  scale_x_continuous(limits = xlim_all) +
  scale_color_viridis_d(name = "Threshold", option = "plasma") +
  labs(x = "Genomic position", y = "HOR Pattern") +
  theme_bw(base_size = 40) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 40),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.ticks.length.y = unit(1, "mm")
  )

# ---------- Save exactly one, depending on --plotv ----------
if (!nzchar(PLOTV) || PLOTV == "V1") {
  out_file <- file.path(OUT_DIR, sprintf("HOR_dotplot_%s.pdf", chr_tag))
  ggsave(out_file, p, height = 10, width = 20, dpi = 300)
  message("[V1] saved: ", out_file)
} else if (PLOTV == "V2") {
  out_file <- file.path(OUT_DIR, sprintf("HOR_dotplot_%s_grdcol.pdf", chr_tag))
  ggsave(out_file, p2, height = 10, width = 20, dpi = 300)
  message("[V2] saved: ", out_file)
} else if (PLOTV == "V3") {
  out_file <- file.path(OUT_DIR, sprintf("HOR_dotplot_%s_pattcor.pdf", chr_tag))
  ggsave(out_file, p3, height = 10, width = 20, dpi = 300)
  message("[V3] saved: ", out_file)
}