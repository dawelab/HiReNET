#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(tidyr)
  library(stringr)
  library(rpart)
  library(data.table)
})

`%!in%` <- Negate(`%in%`)

# ---------------- CLI ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) return(TRUE)
  args[i + 1]
}

INFILE <- get_arg("--input")
OUTDIR <- get_arg("--outdir")
PREFIX <- get_arg("--prefix", "Arabidopsis_CEN178")
BIN    <- as.numeric(get_arg("--bin", "10000"))

if (is.null(INFILE) || is.null(OUTDIR)) {
  stop("Usage: S2_smoothing_bin.R --input FILE --outdir DIR [--prefix NAME] [--bin 10000]\n", call. = FALSE)
}
if (!file.exists(INFILE)) stop("Input file not found: ", INFILE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("[*] INPUT  : ", normalizePath(INFILE))
message("[*] OUTDIR : ", normalizePath(OUTDIR))
message("[*] PREFIX : ", PREFIX)
message("[*] BIN(bp): ", BIN)

# ---------------- Helpers ----------------
# Smooth single-bin spikes: if a bin differs from both neighbors and neighbors agree (and same threshold),
# adopt neighbor class/threshold/prob.
spur_smooth <- function(k, bin_size) {
  bins <- k
  # (re)define end by bin size (don’t trust input end)
  bins$end <- bins$start + bin_size
  bins$lda_pred_smooth <- NA_character_
  bins_smooth <- bins[0, , drop = FALSE]

  for (arr in unique(bins$array)) {
    sub <- bins %>% dplyr::filter(array == arr) %>% arrange(start)
    sub$num <- seq_len(nrow(sub))
    for (j in seq_len(nrow(sub))) {
      if (j > 2 && j < (nrow(sub) - 2)) {
        if (sub$lda_pred[j] != sub$lda_pred[j - 1] &&
            sub$lda_pred[j] != sub$lda_pred[j + 1] &&
            sub$lda_pred[j + 1] == sub$lda_pred[j - 1] &&
            sub$threshold[j + 1] == sub$threshold[j - 1]) {

          sub$lda_pred_smooth[j] <- as.character(sub$lda_pred[j + 1])
          sub$threshold[j]       <- sub$threshold[j + 1]
          sub$lda_pred_pos[j]    <- sub$lda_pred_pos[j + 1]
        }
      }
      if (is.na(sub$lda_pred_smooth[j])) sub$lda_pred_smooth[j] <- as.character(sub$lda_pred[j])
    }
    bins_smooth <- dplyr::bind_rows(bins_smooth, sub)
  }
  bins_smooth
}

# Combine consecutive bins with same class & threshold; merge tiny bins; enforce class priority during overlap cleanup.
comb_over <- function(fil) {
  bins_smooth <- fil
  bins_smooth$tar <- NA_character_
  bins_smooth_comb <- bins_smooth[0, , drop = FALSE]

  for (arr in unique(bins_smooth$array)) {
    sub <- bins_smooth %>% dplyr::filter(array == arr) %>% arrange(start)

    # merge bins with <5 monomers into neighbor
    change <- 1
    while (change > 0) {
      change <- 0
      for (r in seq_len(nrow(sub))) {
        if (sub$num_monomer[r] < 5 && r < nrow(sub) && sub$start[r + 1] != sub$start[r]) {
          sub$start[r + 1]   <- sub$start[r]
          sub$tar[r]         <- "remove"
          sub$num_monomer[r + 1] <- sub$num_monomer[r] + sub$num_monomer[r + 1]
          change <- change + 1
        } else if (sub$num_monomer[r] < 5 && r == nrow(sub) && is.na(sub$tar[r])) {
          sub$end[r - 1]     <- sub$end[r]
          sub$tar[r]         <- "remove"
          sub$num_monomer[r - 1] <- sub$num_monomer[r] + sub$num_monomer[r - 1]
          change <- change + 1
        }
      }
    }

    # keep rows not marked for removal
    sub <- sub[is.na(sub$tar), ]
    sub$num <- seq_len(nrow(sub))

    # group consecutive bins with same smoothed class & threshold
    j <- 2
    tar <- 1
    sub$tar <- NA_character_
    sub$tar[tar] <- "YES"
    while (j <= nrow(sub)) {
      if (sub$lda_pred_smooth[j] == sub$lda_pred_smooth[tar] &&
          sub$start[j] <= sub$end[tar] &&
          sub$end[j]   >= sub$start[tar] &&
          sub$threshold[j] == sub$threshold[tar]) {

        sub$end[tar] <- sub$end[j]
      } else {
        tar <- sub$num[j]
        sub$tar[tar] <- "YES"
      }
      j <- j + 1
    }
    sub2 <- sub[sub$tar %in% "YES", , drop = FALSE]

    # Resolve overlaps giving preference (1) HOR, (2) Order, (3) Disorder
    if (nrow(sub2) > 1) {
      # HOR priority
      for (k in seq_len(nrow(sub2))) {
        if (sub2$lda_pred_smooth[k] == "HOR") {
          if (k == 1) {
            if (nrow(sub2) >= 2) sub2$start[k + 1] <- sub2$end[k]
          } else if (k == nrow(sub2)) {
            sub2$end[k - 1] <- sub2$start[k]
          } else {
            sub2$start[k + 1] <- sub2$end[k]
            sub2$end[k - 1]   <- sub2$start[k]
          }
        }
      }
      # Order priority
      for (k in seq_len(nrow(sub2))) {
        if (sub2$lda_pred_smooth[k] == "Order") {
          if (k == 1) {
            if (nrow(sub2) >= 2) sub2$start[k + 1] <- sub2$end[k]
          } else if (k == nrow(sub2)) {
            sub2$end[k - 1] <- sub2$start[k]
          } else {
            sub2$start[k + 1] <- sub2$end[k]
            sub2$end[k - 1]   <- sub2$start[k]
          }
        }
      }
      # Disorder cleanup
      if (nrow(sub2) >= 2 && sub2$lda_pred_smooth[1] == "Disorder") {
        sub2$end[1] <- sub2$start[2]
      }
      for (k in seq_len(nrow(sub2))) {
        if (sub2$lda_pred_smooth[k] == "Disorder") {
          if (k == 1) {
            if (nrow(sub2) >= 2) sub2$start[k + 1] <- sub2$end[k]
          } else if (k == nrow(sub2)) {
            sub2$end[k - 1] <- sub2$start[k]
          } else {
            sub2$start[k + 1] <- sub2$end[k]
            sub2$end[k - 1]   <- sub2$start[k]
          }
        }
      }

      # regroup in case edges shifted
      j <- 2
      tar <- 1
      sub2$tar <- NA_character_
      sub2 <- sub2[sub2$end > sub2$start, , drop = FALSE]
      sub2$tar[tar] <- "YES"
      sub2$num <- seq_len(nrow(sub2))
      while (j <= nrow(sub2)) {
        if (sub2$lda_pred_smooth[j] == sub2$lda_pred_smooth[tar] &&
            sub2$start[j] <= sub2$end[tar] &&
            sub2$end[j]   >= sub2$start[tar] &&
            sub2$threshold[j] == sub2$threshold[tar]) {

          sub2$end[tar] <- sub2$end[j]
        } else {
          tar <- sub2$num[j]
          sub2$tar[tar] <- "YES"
        }
        j <- j + 1
      }
      sub3 <- sub2[sub2$tar %in% "YES", , drop = FALSE]

      # final simple overlap removal
      if (nrow(sub3) > 1) {
        for (f in seq_len(nrow(sub3) - 1)) {
          if (sub3$end[f] != sub3$start[f + 1]) sub3$start[f + 1] <- sub3$end[f]
        }
      }
    } else {
      sub3 <- sub2
    }

    bins_smooth_comb <- dplyr::bind_rows(bins_smooth_comb, sub3)
  }

  bins_smooth_comb
}

# ---------------- Main ----------------
combined_Multi_bins <- read.table(INFILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# validate expected columns
needed <- c("array","chr","start","end","threshold","num_monomer","dup_monomers",
            "uncon_clust","prop_clust1","prop_clust2","modularity","mean_dist",
            "avg_jacc","max_contract","lda_pred","lda_pred_pos")
missing <- setdiff(needed, colnames(combined_Multi_bins))
if (length(missing)) {
  stop("Missing required columns in input: ", paste(missing, collapse = ", "))
}

combined_Multi_bins2 <- spur_smooth(combined_Multi_bins, BIN)
combined_Multi_bins3 <- comb_over(combined_Multi_bins2)

outfile <- file.path(OUTDIR, paste0(PREFIX, "_fin_bins_combined.txt"))
write.table(combined_Multi_bins3, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

message("[✓] Smoothing/merging complete: ", normalizePath(outfile))