#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rlang)
})

options(scipen = 999)
`%!in%` <- Negate(`%in%`)

# ---------------- CLI ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (i == length(args)) return(TRUE)
  args[i + 1]
}

GROUPS_DIR <- get_arg("--groups-dir")   # dir of *_bins_nam_groups.csv
REVAL_FILE <- get_arg("--reval")        # HOR_re-eval_clusters.txt (has threshold per bin)
COOR_FILE  <- get_arg("--coor")         # all_monomer_bed_inbin.txt (bin monomer coords)
OUTDIR     <- get_arg("--outdir", "horpatterns_out")

if (is.null(GROUPS_DIR) || is.null(REVAL_FILE) || is.null(COOR_FILE)) {
  stop(paste(
    "Usage:",
    "Rscript S4_HORpatterns.R",
    "--groups-dir DIR",
    "--reval HOR_re-eval_clusters.txt",
    "--coor all_monomer_bed_inbin.txt",
    "[--outdir horpatterns_out]"
  ), call. = FALSE)
}

if (!dir.exists(GROUPS_DIR)) stop("groups-dir not found: ", GROUPS_DIR)
if (!file.exists(REVAL_FILE)) stop("reval file not found: ", REVAL_FILE)
if (!file.exists(COOR_FILE))  stop("coor file not found: ", COOR_FILE)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("[*] groups-dir: ", normalizePath(GROUPS_DIR))
message("[*] reval     : ", normalizePath(REVAL_FILE))
message("[*] coor      : ", normalizePath(COOR_FILE))
message("[*] outdir    : ", normalizePath(OUTDIR))

# ---------------- Inputs ----------------
data_files <- dir(GROUPS_DIR, pattern = "_bins_nam_groups\\.csv$", full.names = TRUE)
if (!length(data_files)) stop("No *_bins_nam_groups.csv under: ", GROUPS_DIR)

reval <- read.table(REVAL_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# allow threshold column name to vary (*** fixed braces ***)
if ("threshold" %in% names(reval)) {
  thr_col <- "threshold"
} else if ("clust_val" %in% names(reval)) {
  thr_col <- "clust_val"
} else {
  thr_col <- NA_character_
}
if (is.na(thr_col)) stop("No 'threshold' or 'clust_val' column in reval file.")

# coor: chr array_start array_end pos_start pos_end bin_start bin_end
coor <- read.table(COOR_FILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
if (ncol(coor) < 7) stop("Expected 7 columns in --coor (chr array_start array_end pos_start pos_end bin_start bin_end).")
colnames(coor) <- c("chr","array_start","array_end","pos_start","pos_end","bin_start","bin_end")

coor <- coor %>%
  group_by(chr, array_start, array_end, bin_start, bin_end) %>%
  arrange(chr, array_start, array_end, bin_start, bin_end, pos_start) %>%
  mutate(order = row_number()) %>%
  ungroup()

# ---------------- Config ----------------
clust_nams <- c(LETTERS[1:25], 1:9, letters[1:25])   # labels pool
priv_lab   <- "Z"                                     # private cluster label

# ---------------- Helpers ----------------
pattern_string <- function(dat) {
  # expects columns: nam, clust   (from *_bins_nam_groups.csv)
  clust_uni <- dat %>% dplyr::count(clust, name = "count")
  singles   <- clust_uni$clust[clust_uni$count == 1]
  dat2 <- dat
  if (length(singles) > 0) dat2$clust[dat2$clust %in% singles] <- max(singles)

  ren_tab <- data.frame(clust_val = sort(unique(dat2$clust)))
  ren_tab$letters <- clust_nams[seq_len(nrow(ren_tab))]
  if (length(singles) > 0) ren_tab$letters[ren_tab$clust_val == max(singles)] <- priv_lab

  parts <- stringr::str_split_fixed(dat2$nam, "_", 5)
  dat2$pattern_start <- suppressWarnings(as.numeric(parts[, 4]))
  dat2 <- dat2[order(dat2$pattern_start), , drop = FALSE]
  dat2$letter <- ren_tab$letters[match(dat2$clust, ren_tab$clust_val)]
  dplyr::select(dat2, nam, letter)
}

kmer_counts <- function(seq, kval){
  if (nchar(seq) < kval) return(character(0))
  start <- 1; end <- kval; mers <- character(0)
  while (end <= nchar(seq)) {
    mers <- c(mers, substr(seq, start, end))
    start <- start + 1; end <- end + 1
  }
  mers
}

filtered_kmer_counts <- function(sequence){
  check <- 1; k <- 2
  dat_out <- NULL
  while (check >= 1) {
    out <- kmer_counts(sequence, k)
    k_out <- as.data.frame(table(out), stringsAsFactors = FALSE)
    check <- nrow(k_out[k_out$Freq > 1, ])
    dat_out <- if (k == 2) k_out else if (check > 0) rbind(dat_out, k_out) else dat_out
    k <- k + 1
  }
  if (is.null(dat_out) || !nrow(dat_out)) return(data.frame(out=character(0), Freq=integer(0)))

  # allow overlaps
  for (i in seq_len(nrow(dat_out))) {
    dat_out$Freq[i] <- length(unlist(gregexpr(dat_out$out[i], sequence, fixed = TRUE)))
  }
  # remove 'Z' and singleton kmers
  dat_out_Z <- if (length(grep("Z", dat_out$out)) > 0) dat_out[-grep("Z", dat_out$out), ] else dat_out
  dat_out_Z_sub <- dat_out_Z[dat_out_Z$Freq > 1, , drop = FALSE]

  sum_uniq_lets <- function(x) sum(stringr::str_count(x, LETTERS) > 0)
  max_char_p    <- function(x) max(stringr::str_count(x, LETTERS) / nchar(x))

  dat_out_Z_sub$charuni     <- vapply(as.character(dat_out_Z_sub$out), sum_uniq_lets, numeric(1))
  dat_out_Z_sub$char_prop_m <- vapply(as.character(dat_out_Z_sub$out), max_char_p, numeric(1))
  dat_out_Z_sub$char_num    <- nchar(as.character(dat_out_Z_sub$out))
  dat_out_Z_sub_char <- dat_out_Z_sub[dat_out_Z_sub$charuni > 1 &
                                      dat_out_Z_sub$char_prop_m <= .50 &
                                      dat_out_Z_sub$char_num > 2, , drop = FALSE]

  if (nrow(dat_out_Z_sub_char) > 1) {
    dat_out_Z_sub_char$filt <- NA_character_
    for (i in seq_len(nrow(dat_out_Z_sub_char))) {
      if (is.na(dat_out_Z_sub_char$filt[i])) {
        tar_freq <- dat_out_Z_sub_char$Freq[i]
        tar_num  <- dat_out_Z_sub_char$char_num[i]
        val      <- dat_out_Z_sub_char$out[i]
        dat_out_Z_sub_char <- dat_out_Z_sub_char %>%
          mutate(filt = ifelse(char_num > tar_num & grepl(val, out, fixed = TRUE) & Freq < tar_freq, "filt", filt))
      }
    }
    fin <- dat_out_Z_sub_char[is.na(dat_out_Z_sub_char$filt), , drop = FALSE]
    dat_out_Z_sub2 <- fin[0, ]
    for (i in seq_len(nrow(fin))) {
      tar_freq <- fin$Freq[i]
      tar_num  <- fin$char_num[i]
      val      <- fin$out[i]
      match_tar <- fin %>%
        filter(char_num >= tar_num) %>%
        filter(grepl(val, out, fixed = TRUE)) %>%
        filter(Freq == tar_freq) %>%
        arrange(dplyr::desc(char_num)) %>%
        slice_head(n = 1)
      dat_out_Z_sub2 <- unique(rbind(dat_out_Z_sub2, match_tar))
    }
    return(dat_out_Z_sub2[, c("out","Freq"), drop = FALSE])
  } else {
    return(dat_out_Z_sub_char[, c("out","Freq"), drop = FALSE])
  }
}

bed_convert <- function(bed_fil, patt){
  # bed_fil must have a character column in the first position
  if (!nrow(bed_fil)) return(bed_fil)
  v <- as.character(bed_fil[[1]])
  as_bed <- data.frame(V1=character(0), V2=integer(0), V3=integer(0), stringsAsFactors = FALSE)
  for (i in seq_along(v)) {
    token  <- as.character(v[i])
    coords <- unlist(gregexpr(token, patt, fixed = TRUE))
    coords <- coords[coords > 0]
    if (!length(coords)) next
    co_sub <- data.frame(
      V1 = token,
      V2 = coords,
      V3 = coords + nchar(token) - 1,
      stringsAsFactors = FALSE
    )
    as_bed <- rbind(as_bed, co_sub)
  }
  as_bed
}

HOR_patt_check <- function(bed_fil, patt){
  as_bed <- bed_convert(bed_fil, patt)
  if (!nrow(as_bed)) return(as_bed)
  as_bed_or <- as_bed[order(as_bed$V2), ]
  changes <- 1
  while (changes > 0) {
    changes <- 0
    end <- if (nrow(as_bed_or) <= 2) 2 else nrow(as_bed_or) - 1
    for (i in 2:end) {
      if (as_bed_or[i, 2] <= as_bed_or[i - 1, 3]) {
        as_bed_or[i, 2] <- as_bed_or[i - 1, 3] + 1
        changes <- changes + 1
      }
    }
    as_bed_or <- as_bed_or[as_bed_or$V3 - as_bed_or$V2 > 0, ]
  }
  as_bed_or$chr_num    <- nchar(as_bed_or$V1)
  as_bed_or$patt_coord <- as_bed_or$V3 - as_bed_or$V2 + 1
  as_bed_or2 <- as_bed_or %>%
    filter(patt_coord >= chr_num * .5) %>%
    group_by(V1) %>% filter(n() >= 1) %>% dplyr::count(V1, name = "n")
  if (!nrow(as_bed_or2)) return(as_bed_or2)
  as_bed_2 <- bed_convert(as_bed_or2["V1"], patt)
  as_bed_2[order(as_bed_2$V2), , drop = FALSE]
}

purity_calc <- function(as_bed_2_or, patt){
  if (!nrow(as_bed_2_or)) return(0)
  as_bed_2_or$V2 <- as.numeric(as_bed_2_or$V2)
  as_bed_2_or$V3 <- as.numeric(as_bed_2_or$V3)
  as_bed_int <- data.frame(V1 = as_bed_2_or$V2[1], V2 = as_bed_2_or$V3[1])
  start <- 1
  for (i in 2:nrow(as_bed_2_or)) {
    if (as_bed_2_or[i, 2] <= as_bed_int[start, 2]) {
      as_bed_int[start, 2] <- as_bed_2_or[i, 3]
    } else {
      start <- start + 1
      as_bed_int[start, 1] <- as_bed_2_or[i, 2]
      as_bed_int[start, 2] <- as_bed_2_or[i, 3]
    }
  }
  as_bed_int$V3 <- as_bed_int$V2 - as_bed_int$V1 + 1
  len <- sum(as_bed_int$V3)
  len / nchar(patt)
}

group_letter_monomer <- function(pattern_out, HOR_bed_filt,
                                 outdir_root = ".",
                                 only_keep_if_letter_in_pattern = TRUE) {
  dir.create(outdir_root, showWarnings = FALSE, recursive = TRUE)
  patterns_map <- HOR_bed_filt %>%
    select(chr, array_start, array_end, bin_start, bin_end, threshold, patterns) %>%
    distinct()

  join_keys <- intersect(
    c("chr","array_start","array_end","bin_start","bin_end","threshold"),
    intersect(names(pattern_out), names(patterns_map))
  )

  df <- pattern_out %>% left_join(patterns_map, by = join_keys)

  if (only_keep_if_letter_in_pattern && "patterns" %in% names(df)) {
    df <- df %>% filter(mapply(function(l, p) grepl(l, p, fixed = TRUE), letter, patterns))
  }

  out_tbl <- df %>%
    mutate(
      th_str  = sprintf("%.2f", as.numeric(threshold)),
      fa_name = paste0(nam, ".fa"),
      subdir  = paste(chr, array_start, array_end, bin_start, bin_end, sep = "_"),
      filename = paste(chr, array_start, array_end, bin_start, bin_end, letter, th_str, sep = "_"),
      filename = paste0(filename, ".list.txt"),
      outdir   = file.path(outdir_root, subdir),
      outpath  = file.path(outdir, filename)
    ) %>%
    select(outdir, outpath, fa_name) %>%
    distinct()

  lst <- split(out_tbl$fa_name, out_tbl$outpath)
  purrr::iwalk(lst, ~ {
    dir.create(dirname(.y), showWarnings = FALSE, recursive = TRUE)
    writeLines(.x, .y)
  })
  invisible(out_tbl)
}

# ---------------- Outputs init ----------------
string_out <- as.data.frame(matrix(nrow = length(data_files), ncol = 9))
colnames(string_out) <- c("chr","array_start","array_end","bin_start","bin_end","threshold","num_monomers","pattern","purity")

HOR_bed_filt_all <- as.data.frame(matrix(nrow = 0, ncol = 11))
colnames(HOR_bed_filt_all) <- c("chr","array_start","array_end","bin_start","bin_end","threshold",
                                "num_monomers","pattern_start","pattern_end","patterns","purity")

# ---------------- Main loop ----------------
for (j in seq_along(data_files)) {
  dat <- read.csv(data_files[j], stringsAsFactors = FALSE)
  bx  <- basename(data_files[j])

  # Accept either "... .blat.sub_bins_nam_groups.csv" or "... _bins_nam_groups.csv"
  m <- stringr::str_match(
    bx,
    "^(chr[^_]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)(?:\\.blat\\.sub)?_bins_nam_groups\\.csv$"
  )
  if (any(is.na(m))) stop("Unexpected groups filename: ", bx)

  string_out$chr[j]         <- m[2]
  string_out$array_start[j] <- as.numeric(m[3])
  string_out$array_end[j]   <- as.numeric(m[4])
  string_out$bin_start[j]   <- as.numeric(m[5])
  string_out$bin_end[j]     <- as.numeric(m[6])

  array_key <- paste0(string_out$array_start[j], "-", string_out$array_end[j])

  idx <- reval$chr == string_out$chr[j] &
         reval$start == string_out$bin_start[j] &
         reval$end   == string_out$bin_end[j] &
         (!"array" %in% names(reval) | reval$array == array_key)

  string_out$threshold[j]    <- if (any(idx)) reval[idx, thr_col][1] else NA_real_
  string_out$num_monomers[j] <- nrow(dat)

  # pattern letters by monomer order
  pattern_out <- pattern_string(dat) %>%
    separate(
      nam,
      into = c("chr", "array_start", "array_end", "monomer_start", "monomer_end"),
      sep = "_", remove = FALSE, fill = "right"
    ) %>%
    mutate(across(c(array_start, array_end, monomer_start, monomer_end),
                  ~ suppressWarnings(as.numeric(.))),
           bin_start = string_out$bin_start[j],
           bin_end   = string_out$bin_end[j],
           threshold = string_out$threshold[j])

  string_out$pattern[j] <- paste(pattern_out$letter, collapse = "")

  # Find HOR patterns
  HOR_patterns <- filtered_kmer_counts(string_out$pattern[j])
  if (!is.data.frame(HOR_patterns)) {
    HOR_patterns <- as.data.frame(HOR_patterns, stringsAsFactors = FALSE)
  }
  names(HOR_patterns) <- c("out","Freq")[seq_len(ncol(HOR_patterns))]

  if (nrow(HOR_patterns) > 0) {
    HOR_bed_filt <- HOR_patt_check(HOR_patterns["out"], string_out$pattern[j])
  } else {
    HOR_bed_filt <- HOR_patterns[0, ]
  }

  if (nrow(HOR_bed_filt) > 0) {
    string_out$purity[j] <- purity_calc(HOR_bed_filt, string_out$pattern[j])

    colnames(HOR_bed_filt) <- c("patterns","pattern_start","pattern_end")
    HOR_bed_filt$chr          <- string_out$chr[j]
    HOR_bed_filt$array_start  <- string_out$array_start[j]
    HOR_bed_filt$array_end    <- string_out$array_end[j]
    HOR_bed_filt$bin_start    <- string_out$bin_start[j]
    HOR_bed_filt$bin_end      <- string_out$bin_end[j]
    HOR_bed_filt$threshold    <- string_out$threshold[j]
    HOR_bed_filt$num_monomers <- string_out$num_monomers[j]
    HOR_bed_filt$purity       <- string_out$purity[j]

    HOR_bed_filt <- HOR_bed_filt %>%
      select(chr,array_start,array_end,bin_start,bin_end,threshold,
             num_monomers,pattern_start,pattern_end,patterns,purity)

    HOR_bed_filt_all <- rbind(HOR_bed_filt_all, HOR_bed_filt)
    HOR_bed_filt_all$pattern_end <- HOR_bed_filt_all$pattern_start +
                                    nchar(HOR_bed_filt_all$patterns) - 1

    # write nam↔letter map per bin
    outdir1 <- file.path(OUTDIR, "mergebin_string_outputs")
    dir.create(outdir1, showWarnings = FALSE, recursive = TRUE)
    outfile <- paste(string_out$chr[j], string_out$array_start[j], string_out$array_end[j],
                     string_out$bin_start[j], string_out$bin_end[j], sep = "_")
    outpath <- file.path(outdir1, paste0(outfile, ".out"))
    write.table(pattern_out[, c("nam","letter")], outpath, col.names = FALSE, row.names = FALSE, quote = FALSE)

    # write group lists by letter
    outdir2 <- file.path(OUTDIR, "group_HOR_monomers")
    dir.create(outdir2, showWarnings = FALSE, recursive = TRUE)
    invisible(group_letter_monomer(pattern_out, HOR_bed_filt, outdir_root = outdir2))
  } else {
    string_out$purity[j] <- 0
  }
}

# ---------------- Save tables ----------------
out_string <- file.path(OUTDIR, "string_out.txt")
out_horbed <- file.path(OUTDIR, "HOR_bed.txt")
write.table(HOR_bed_filt_all, out_horbed, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(string_out,      out_string, sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------- Map to coordinates ----------------
keys <- c("chr","array_start","array_end","bin_start","bin_end")

coor_start <- coor %>%
  transmute(
    !!!syms(keys),
    pattern_start = order,
    Monomer_first_start = pos_start,
    Monomer_first_end   = pos_end
  )

coor_end <- coor %>%
  transmute(
    !!!syms(keys),
    pattern_end = order,
    Monomer_end_start = pos_start,
    Monomer_end_end   = pos_end
  )

HOR_bed <- HOR_bed_filt_all

out <- HOR_bed %>%
  select(all_of(c(keys, "threshold","num_monomers","pattern_start","pattern_end","patterns","purity"))) %>%
  left_join(coor_start, by = c(keys, "pattern_start")) %>%
  left_join(coor_end,   by = c(keys, "pattern_end"))

out_horcoor <- file.path(OUTDIR, "HOR_coor_bed.txt")
write.table(out, out_horcoor, sep = "\t", quote = FALSE, row.names = FALSE)

message("[✓] Done.")
message("  - ", normalizePath(out_string))
message("  - ", normalizePath(out_horbed))
message("  - ", normalizePath(out_horcoor))
message("  - per-bin maps: ", normalizePath(file.path(OUTDIR, "mergebin_string_outputs")))
message("  - group lists:  ", normalizePath(file.path(OUTDIR, "group_HOR_monomers")))