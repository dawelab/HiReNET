#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ---------------- CLI parsing ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) return(TRUE)
  args[i + 1]
}

HOR_FILE      <- get_arg("--hor")
OUTDIR        <- get_arg("--outdir", ".")
BINWIDTH_HOR1 <- as.numeric(get_arg("--binwidth-hor", "170"))   # for the first HOR histogram
BINWIDTH_HOR2 <- as.numeric(get_arg("--binwidth-hor2", "10"))   # for the “spiky” zoomed plot
XMAX_HOR2     <- as.numeric(get_arg("--xmax-hor", "2000"))

if (is.null(HOR_FILE)) {
  stop(paste(
    "Usage:",
    "Rscript S5_plot_cli.R --hor HOR_coor_bed.txt [--outdir plots_out] [--binwidth-hor 170] [--binwidth-hor2 10] [--xmax-hor 2000]"
  ), call. = FALSE)
}

if (!file.exists(HOR_FILE)) stop("HOR file not found: ", HOR_FILE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("[*] HOR file: ", normalizePath(HOR_FILE))
message("[*] outdir  : ", normalizePath(OUTDIR))

# ---------------- Load & prep ----------------
# Expected columns (from your pipeline):
# chr array_start array_end bin_start bin_end threshold num_monomers
# pattern_start pattern_end patterns purity Monomer_first_start Monomer_first_end Monomer_end_start Monomer_end_end
hor <- suppressMessages(read_tsv(HOR_FILE, show_col_types = FALSE))

# Robustness: ensure required columns exist
needed <- c("patterns", "Monomer_first_start", "Monomer_end_end")
missing <- setdiff(needed, names(hor))
if (length(missing)) {
  stop("Missing columns in HOR file: ", paste(missing, collapse = ", "))
}

hor <- hor %>%
  mutate(
    monomer_size   = nchar(patterns),
    HOR_length     = abs(as.numeric(Monomer_end_end) - as.numeric(Monomer_first_start)),
    HOR_length_est = monomer_size * 200
  )

# Filter (like your code): keep HORs shorter than their rough estimate
hor_fil <- hor %>% filter(HOR_length < HOR_length_est)

# ---------------- Plots ----------------
# 1) Monomer size frequency
Mon_fre <- ggplot(hor_fil, aes(x = monomer_size)) +
  geom_bar(fill = "steelblue") +
  scale_x_continuous(breaks = seq(0, max(hor_fil$monomer_size, na.rm = TRUE), by = 3)) +
  labs(
    x = "Monomer size",
    y = "Frequency",
    title = "HOR monomer size frequency"
  ) +
  theme_classic(base_size = 14)

# 2) HOR length frequency (binned, broad)
Len_fre1 <- ggplot(hor_fil, aes(x = HOR_length)) +
  geom_histogram(binwidth = BINWIDTH_HOR1, fill = "black", color = "white") +
  scale_x_continuous(
    breaks = seq(0, max(hor$HOR_length, na.rm = TRUE), by = 1000)
  ) +
  labs(
    x = "HOR length (bp)",
    y = "Frequency",
    title = "HOR length frequency"
  ) +
  theme_classic(base_size = 14)

# 3) HOR length frequency (spiky, zoomed to [0, XMAX_HOR2])
Len_fre2 <- ggplot(hor, aes(x = HOR_length)) +
  geom_histogram(
    binwidth = BINWIDTH_HOR2,
    color    = "black",
    fill     = NA,
    linewidth = 1,
    boundary = 0,
    closed   = "left"
  ) +
  coord_cartesian(xlim = c(0, XMAX_HOR2)) +
  scale_x_continuous(breaks = c(0, XMAX_HOR2/2, XMAX_HOR2)) +
  labs(
    x = "HOR length (bp)",
    y = "Frequency",
    title = "HOR length frequency (zoomed)"
  ) +
  theme_classic(base_size = 14)

# ---------------- Save ----------------
pdf1 <- file.path(OUTDIR, "Monomersize_freq.pdf")
pdf2 <- file.path(OUTDIR, "HOR_length_freq.pdf")
pdf3 <- file.path(OUTDIR, "HOR_length_freq_zoomed.pdf")

ggsave(pdf1, Mon_fre, height = 8, width = 10, dpi = 300)
ggsave(pdf2, Len_fre1, height = 8, width = 10, dpi = 300)
ggsave(pdf3, Len_fre2, height = 8, width = 10, dpi = 300)

message("[✓] Wrote: ", normalizePath(pdf1))
message("[✓] Wrote: ", normalizePath(pdf2))
message("[✓] Wrote: ", normalizePath(pdf3))