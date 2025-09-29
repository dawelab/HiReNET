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

# ---------- read data (assign thresholds by order 91–99) ----------
data_files <- dir(STRING_DIR, pattern = "_HOR_patterns_bed\\.csv$", full.names = TRUE)
if (!length(data_files)) {
  message("No *_HOR_patterns_bed.csv found under ", STRING_DIR)
  quit(save = "no", status = 0)
}

data_files <- sort(data_files)  # stable order
thresholds <- 91:99
if (length(data_files) != length(thresholds)) {
  stop(sprintf("Expected %d files for thresholds 91–99, but found %d.\nFiles:\n%s",
               length(thresholds), length(data_files),
               paste(basename(data_files), collapse = "\n")))
}

normalize_types <- function(d) {
  d %>%
    mutate(
      chr          = as.character(chr),
      array_start  = as.numeric(array_start),
      array_end    = as.numeric(array_end),
      bin_start    = as.numeric(bin_start),
      bin_end      = as.numeric(bin_end),
      pattern_start= as.numeric(pattern_start),
      pattern_end  = as.numeric(pattern_end)
    )
}

dat_list <- lapply(seq_along(data_files), function(i) {
  d <- read.csv(data_files[i], stringsAsFactors = FALSE)
  if (!nrow(d)) return(dplyr::tibble())
  normalize_types(d) %>% mutate(threshold = thresholds[i])
})

dat_list <- Filter(function(x) nrow(x) > 0, dat_list)
if (!length(dat_list)) {
  message("All *_HOR_patterns_bed.csv are empty. Nothing to plot.")
  quit(save = "no", status = 0)
}

# -------- combine & derive fields ----------
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
    pos = array_start + (pattern_start + pattern_end) / 2
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

# -------- dot plot, no filtering ----------
p <- ggplot(df_all, aes(x = pos, y = order_id, color = threshold)) +
  geom_point(alpha = 0.85, size = 1) +
  scale_y_continuous(
    breaks = as.numeric(names(lab_map)), labels = lab_map,
    limits = c(0, max(df_all$order_id, na.rm = TRUE)),
    expand = expansion(mult = c(0, 0), add = 0)
  ) +
  scale_x_continuous(limits = xlim_all) +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "Genomic position", y = "HOR Pattern", color = "Threshold") +
  theme_bw(base_size = 40) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 40))

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

ggsave(file.path(OUT_DIR, "shared_HOR.pdf"),             p,        height = 40, width = 40, dpi = 300)
ggsave(file.path(OUT_DIR, "shared_HOR_mon_size.pdf"),    p_size,   height = 10, width = 15, dpi = 300)
ggsave(file.path(OUT_DIR, "shared_HOR_counts_10kb.pdf"), p_count,  height = 10, width = 15, dpi = 300)

# -------- dot plot (threshold > 95) ----------
df_all_fil <- df_all %>% filter(threshold > 95)
if (nrow(df_all_fil)) {
  df_all_fil <- df_all_fil %>%
    arrange(threshold, HOR_pattern) %>%
    group_by(threshold, HOR_pattern) %>%
    mutate(order_id = cur_group_id()) %>%
    ungroup()

  lab_map_fil <- setNames(df_all_fil$HOR_pattern, df_all_fil$order_id)
  xlim_fil <- c(min(df_all_fil$array_start, na.rm = TRUE),
                max(df_all_fil$array_end,   na.rm = TRUE))

  p_fil <- ggplot(df_all_fil, aes(x = pos, y = order_id, color = threshold)) +
    geom_point(alpha = 0.85, size = 1) +
    scale_y_continuous(
      breaks = as.numeric(names(lab_map_fil)), labels = lab_map_fil,
      limits = c(0, max(df_all_fil$order_id, na.rm = TRUE)),
      expand = expansion(mult = c(0, 0), add = 0)
    ) +
    scale_x_continuous(limits = xlim_fil) +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "Genomic position", y = "HOR Pattern", color = "Threshold") +
    theme_bw(base_size = 25) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.y = element_text(size = 25),
          legend.text  = element_text(size = 15))

  p_fil2 <- ggplot(df_all_fil, aes(x = pos, y = order_id, color = threshold)) +
    geom_point(alpha = 0.85, size = 1) +
    scale_y_continuous(
      limits = c(0, max(df_all_fil$order_id, na.rm = TRUE)),
      expand = expansion(mult = c(0, 0), add = 0)
    ) +
    scale_x_continuous(limits = xlim_fil) +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "Genomic position", y = "HOR Pattern", color = "Threshold") +
    theme_bw(base_size = 25) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.text  = element_text(size = 15))

  p_fil_size <- ggplot(df_all_fil, aes(x = mon_size)) +
    geom_bar() +
    scale_x_continuous(breaks = c(3, 6, 9, 12, 15)) +
    labs(x = "Number of monomers in HOR pattern ", y = "Frequency") +
    theme_bw(base_size = 25) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggsave(file.path(OUT_DIR, "shared_HOR_fil_thr95.pdf"),          p_fil,      height = 15, width = 20, dpi = 300)
  ggsave(file.path(OUT_DIR, "shared_HOR_fil_thr95_id.pdf"),       p_fil2,     height = 15, width = 20, dpi = 300)
  ggsave(file.path(OUT_DIR, "shared_HOR_mon_size_fil_thr95.pdf"), p_fil_size, height = 10, width = 15, dpi = 300)
} else {
  message("No rows with threshold > 95; skipping *_thr95 plots.")
}

# -------- dot plot (count > 5) ----------
df_all_fil2 <- df_all %>% filter(count > 5)
if (nrow(df_all_fil2)) {
  df_all_fil2 <- df_all_fil2 %>%
    arrange(threshold, HOR_pattern) %>%
    group_by(threshold, HOR_pattern) %>%
    mutate(order_id = cur_group_id()) %>%
    ungroup()

  lab_map_fil2 <- setNames(df_all_fil2$HOR_pattern, df_all_fil2$order_id)
  xlim_fil2 <- c(min(df_all_fil2$array_start, na.rm = TRUE),
                 max(df_all_fil2$array_end,   na.rm = TRUE))

  p_filc <- ggplot(df_all_fil2, aes(x = pos, y = order_id, color = threshold)) +
    geom_point(alpha = 0.85, size = 1) +
    scale_y_continuous(
      breaks = as.numeric(names(lab_map_fil2)), labels = lab_map_fil2,
      limits = c(0, max(df_all_fil2$order_id, na.rm = TRUE)),
      expand = expansion(mult = c(0, 0), add = 0)
    ) +
    scale_x_continuous(limits = xlim_fil2) +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "Genomic position", y = "HOR Pattern", color = "Threshold") +
    theme_bw(base_size = 25) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.title.y = element_text(size = 25),
          legend.text  = element_text(size = 15))

  p_fils <- ggplot(df_all_fil2, aes(x = mon_size)) +
    geom_bar() +
    scale_x_continuous(breaks = c(3, 6, 9, 12, 15)) +
    labs(x = "Number of monomers in HOR pattern ", y = "Frequency") +
    theme_bw(base_size = 25) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggsave(file.path(OUT_DIR, "shared_HOR_fil_count5.pdf"),          p_filc, height = 15, width = 20, dpi = 300)
  ggsave(file.path(OUT_DIR, "shared_HOR_mon_size_fil_count5.pdf"), p_fils, height = 10, width = 15, dpi = 300)
} else {
  message("No rows with count > 5; skipping *count5 plots.")
}