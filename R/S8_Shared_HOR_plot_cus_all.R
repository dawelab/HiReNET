#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
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

STRING_DIR <- file.path(OUT_DIR, "string_share_HOR")
if (!dir.exists(STRING_DIR)) stop("Not found: ", STRING_DIR)

# ---------- read data (RELAXED: infer thresholds, skip empties) ----------
data_files <- dir(STRING_DIR, pattern = "_HOR_patterns_bed\\.csv$", full.names = TRUE)
if (!length(data_files)) {
  message("No *_HOR_patterns_bed.csv found under ", STRING_DIR)
  quit(save = "no", status = 0)
}
data_files <- sort(data_files)

# try to infer threshold (e.g., "…_93_bins_nam_groups_HOR_patterns_bed.csv")
extract_thr <- function(p) {
  b <- basename(p)
  m <- regmatches(b, regexec("_(\\d+)_bins_nam_groups_HOR_patterns_bed\\.csv$", b))[[1]]
  if (length(m) >= 2) as.integer(m[2]) else NA_integer_
}
thr_from_name <- vapply(data_files, extract_thr, integer(1))
if (any(is.na(thr_from_name))) {
  # fallback: assign sequentially from 90 in file order
  thr_from_name <- seq.int(90, length.out = length(data_files))
}

normalize_types <- function(d) {
  d %>%
    mutate(
      chr           = as.character(chr),
      array_start   = as.numeric(array_start),
      array_end     = as.numeric(array_end),
      bin_start     = as.numeric(bin_start),
      bin_end       = as.numeric(bin_end),
      pattern_start = as.numeric(pattern_start),
      pattern_end   = as.numeric(pattern_end)
    )
}

dat_list <- lapply(seq_along(data_files), function(i) {
  d <- tryCatch(read.csv(data_files[i], stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)                 # <-- skip empty
  d <- normalize_types(d)
  if (!nrow(d)) return(NULL)                               # <-- skip if all-NA after coerce
  d$threshold <- thr_from_name[i]
  d
})
dat_list <- Filter(Negate(is.null), dat_list)
if (!length(dat_list)) {
  message("All *_HOR_patterns_bed.csv are empty. Nothing to plot.")
  quit(save = "no", status = 0)
}

# -------- combine & derive fields (unchanged) ----------
df <- dplyr::bind_rows(dat_list)

df_fre <- df %>%
  group_by(HOR_pattern, threshold) %>%
  summarise(count = n(), .groups = "drop")

df_all <- left_join(df, df_fre, by = c("HOR_pattern","threshold")) %>%
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
    pos = (pattern_start + pattern_end) / 2
  )

if (!nrow(df_all)) {
  message("No rows to plot after combining inputs. Exiting.")
  quit(save = "no", status = 0)
}

df_all$HOR_pattern <- str_replace_all(df_all$HOR_pattern, "\\|", " ")
df_all <- df_all %>% mutate(mon_size = stringr::str_count(HOR_pattern, "\\S+"))

lab_map <- setNames(df_all$HOR_pattern, df_all$order_id)

# -------- dot plot, no filtering (unchanged) ----------
p <- ggplot(df_all, aes(x = pos, y = order_id, color = threshold)) +
  geom_point(alpha = 0.85, size = 2) +
  facet_wrap(~chr, scales = "free_x", ncol = 1,strip.position = "right") +
  scale_y_continuous(
    breaks = as.numeric(names(lab_map)), labels = lab_map,
    limits = c(0, max(df_all$order_id, na.rm = TRUE)),
    expand = expansion(mult = c(0, 0), add = 0)
  ) +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Genomic position", y = "HOR Pattern", color = "Threshold") +
  theme_bw(base_size = 30) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 30)
  )

p_size <- ggplot(df_all, aes(x = mon_size)) +
  geom_bar() +
  scale_x_continuous(breaks = c(3, 6, 9, 12, 15)) +
  labs(x = "Number of monomers in shared HOR pattern ", y = "Frequency") +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p_count <- ggplot(df_all, aes(x = pos)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(
    breaks = seq(min(df_all$pos, na.rm = TRUE),
                 max(df_all$pos, na.rm = TRUE), by = 500000)
  ) +
  labs(x = "Genomic position", y = "Shared HOR patterns counts per 10kb window") +
  theme_bw(base_size = 25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "shared_HOR.pdf"),             p,        height = 50, width = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "shared_HOR_mon_size.pdf"),    p_size,   height = 10, width = 15, dpi = 300)
ggsave(file.path(OUT_DIR, "shared_HOR_counts_10kb.pdf"), p_count,  height = 10, width = 15, dpi = 300)

