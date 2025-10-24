#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
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
OUT_DIR     <- argval("--outdir",  default = NULL)
LETTER_DIR  <- argval("--letters", default = NULL)
if (is.null(OUT_DIR)    || OUT_DIR    == "") stop("Missing --outdir")
if (is.null(LETTER_DIR) || LETTER_DIR == "") stop("Missing --letters")

OUT_DIR    <- normalizePath(OUT_DIR,    mustWork = TRUE)
LETTER_DIR <- normalizePath(LETTER_DIR, mustWork = TRUE)

CSV_DIR    <- file.path(OUT_DIR, "csv_share_HOR")        # from S1
STRING_DIR <- file.path(OUT_DIR, "string_share_HOR")     # outputs here
if (!dir.exists(CSV_DIR)) {
  stop("Expected CSV directory from S1 not found: ", CSV_DIR)
}
if (!dir.exists(STRING_DIR)) dir.create(STRING_DIR, recursive = TRUE, showWarnings = FALSE)

message("[*] CSV input dir : ", CSV_DIR)
message("[*] Letter dir    : ", LETTER_DIR)
message("[*] Output dir    : ", STRING_DIR)

# ---------- labeling config (unchanged) ----------
base_alphabet <- c(LETTERS[1:25], letters[1:25])  # A–Y, a–y
clust_nams    <- base_alphabet
priv_lab      <- "Z"
fallback_labels <- function(n) paste0("C", seq_len(n))  # C1, C2, ...

# ---------- helpers ----------
# parse 2-col letter files: "<chr>_<array_start>_<array_end>_<mon_start>_<mon_end> <LETTER>"
read_letters_2col <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  df <- tryCatch(read.table(fp, header = FALSE, sep = "", stringsAsFactors = FALSE,
                            quote = "", comment.char = "", fill = TRUE),
                 error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2L || !nrow(df)) return(NULL)
  colnames(df)[1:2] <- c("id", "letter")
  parts <- strsplit(df$id, "_", fixed = TRUE)
  ok <- lengths(parts) == 5L
  if (!any(ok)) return(NULL)
  parts <- parts[ok]; df <- df[ok, , drop = FALSE]
  mat <- do.call(rbind, parts)
  colnames(mat) <- c("chr","array_start","array_end","monomer_start","monomer_end")
  out <- data.frame(
    monomer       = df$id,
    chr           = mat[, "chr"],
    array_start   = mat[, "array_start"],
    array_end     = mat[, "array_end"],
    monomer_start = as.numeric(mat[, "monomer_start"]),
    monomer_end   = as.numeric(mat[, "monomer_end"]),
    letter        = df$letter,
    stringsAsFactors = FALSE
  )
  out
}

pattern_string_v2 <- function(dat) {
  parts <- stringr::str_split_fixed(dat$nam, "_", 7)
  new_dat <- data.frame(
    chr         = parts[, 1],
    array_start = parts[, 2],
    array_end   = parts[, 3],
    bin_start   = parts[, 4],
    bin_end     = parts[, 5],
    letter      = parts[, 6],
    thr         = parts[, 7],
    clust       = dat$clust,
    stringsAsFactors = FALSE
  )
  new_dat$id <- paste(new_dat$chr, new_dat$array_start, new_dat$array_end,
                      new_dat$bin_start, new_dat$bin_end, sep = "_")

  clust_uni <- new_dat |>
    dplyr::group_by(clust) |>
    dplyr::summarise(count = dplyr::n_distinct(id), .groups = "drop")

  singles  <- clust_uni$clust[clust_uni$count == 1]
  non_sing <- sort(setdiff(unique(new_dat$clust), singles))

  if (length(non_sing) <= length(clust_nams)) {
    label_pool <- clust_nams[seq_along(non_sing)]
  } else {
    label_pool <- c(clust_nams, fallback_labels(length(non_sing) - length(clust_nams)))
  }

  map_non <- tibble::tibble(clust = non_sing, label_non = label_pool)

  new_dat |>
    dplyr::left_join(map_non, by = "clust") |>
    dplyr::mutate(
      new_letter = dplyr::case_when(
        clust %in% singles ~ priv_lab,
        TRUE               ~ label_non
      )
    ) |>
    dplyr::select(-label_non)
}

.tokenize <- function(x) {
  if (length(x) == 1L) strsplit(trimws(x), "\\s+")[[1]] else as.character(x)
}
kmer_counts <- function(seq, kval) {
  tokens <- .tokenize(seq); n <- length(tokens)
  if (n < kval) return(character(0))
  sapply(seq_len(n - kval + 1L), function(i) paste(tokens[i:(i + kval - 1L)], collapse = "|"))
}
.count_nonoverlap_tokens <- function(needle_tokens, hay_tokens) {
  m <- length(needle_tokens); n <- length(hay_tokens)
  if (n < m) return(0L)
  i <- 1L; cnt <- 0L
  while (i <= n - m + 1L) {
    if (all(hay_tokens[i:(i + m - 1L)] == needle_tokens)) { cnt <- cnt + 1L; i <- i + m }
    else i <- i + 1L
  }
  cnt
}
.contains_subseq <- function(long_tokens, short_tokens) {
  m <- length(short_tokens); n <- length(long_tokens)
  if (m == 0L || n < m) return(FALSE)
  any(sapply(seq_len(n - m + 1L),
             function(i) all(long_tokens[i:(i + m - 1L)] == short_tokens)))
}
filtered_kmer_counts <- function(sequence) {
  tokens <- .tokenize(sequence)
  if (!length(tokens)) {
    return(data.frame(out=character(), Freq=integer(),
                      charuni=integer(), char_prop_m=numeric(), char_num=integer(),
                      stringsAsFactors = FALSE))
  }
  check <- 1L; k <- 2L; dat_out <- NULL
  while (check >= 1L) {
    out <- kmer_counts(tokens, k); if (!length(out)) break
    k_out <- as.data.frame(table(out), stringsAsFactors = FALSE)
    names(k_out) <- c("out", "Freq")
    check <- sum(k_out$Freq > 1L)
    if (k == 2L) dat_out <- k_out else if (check > 0) dat_out <- rbind(dat_out, k_out) else break
    k <- k + 1L
  }
  if (is.null(dat_out) || nrow(dat_out) == 0L) {
    return(data.frame(out=character(), Freq=integer(),
                      charuni=integer(), char_prop_m=numeric(), char_num=integer(),
                      stringsAsFactors = FALSE))
  }
  split_pat_all <- strsplit(dat_out$out, "\\|"); tokens_all <- .tokenize(sequence)
  dat_out$Freq <- vapply(split_pat_all,
                         function(pt) .count_nonoverlap_tokens(pt, tokens_all), integer(1L))
  has_Z   <- vapply(split_pat_all, function(v) any(v == "Z"), logical(1))
  dat_out <- dat_out[dat_out$Freq > 1L & !has_Z, , drop = FALSE]
  if (nrow(dat_out) == 0L) {
    return(data.frame(out=character(), Freq=integer(),
                      charuni=integer(), char_prop_m=numeric(), char_num=integer(),
                      stringsAsFactors = FALSE))
  }
  split_pat <- strsplit(dat_out$out, "\\|")
  dat_out$charuni  <- vapply(split_pat, function(v) length(unique(v)), integer(1L))
  dat_out$char_num <- vapply(split_pat, length, integer(1L))
  dat_out$char_prop_m <- vapply(split_pat, function(v) { tv <- table(v); max(as.numeric(tv)) / length(v) }, numeric(1L))
  dat_out <- subset(dat_out, charuni > 1L & char_prop_m <= 0.50 & char_num > 2L)
  rownames(dat_out) <- NULL
  dat_out
}
filtered_kmer_counts_2 <- function(df) {
  if (is.null(df) || nrow(df) <= 1L) return(df)
  df$filt <- NA_character_; split_out <- strsplit(df$out, "\\|")
  for (i in seq_len(nrow(df))) {
    if (!is.na(df$filt[i]) && df$filt[i] == "filt") next
    tar_freq <- df$Freq[i]; tar_num <- df$char_num[i]; needle <- split_out[[i]]
    for (j in seq_len(nrow(df))) {
      if (df$char_num[j] > tar_num &&
          .contains_subseq(split_out[[j]], needle) &&
          df$Freq[j] < tar_freq) df$filt[j] <- "filt"
    }
  }
  fin_idx <- which(is.na(df$filt) | df$filt != "filt")
  fin <- df[fin_idx, , drop = FALSE]; split_fin <- split_out[fin_idx]
  keep_rows <- logical(nrow(fin))
  for (i in seq_len(nrow(fin))) {
    tar_freq <- fin$Freq[i]; needle <- split_fin[[i]]
    candidates <- which(fin$Freq == tar_freq &
                          vapply(split_fin, .contains_subseq, logical(1), short_tokens = needle))
    best <- candidates[which.max(fin$char_num[candidates])]; keep_rows[best] <- TRUE
  }
  unique(fin[keep_rows, , drop = FALSE])
}
bed_convert <- function(bed_fil, patt, starts, ends) {
  tokens <- .tokenize(patt)
  as_bed <- data.frame(V1=character(), V2=integer(), V3=integer(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(bed_fil))) {
    pat_str <- as.character(bed_fil[i, 1])
    pat_tokens <- strsplit(trimws(pat_str), "\\s+|\\|")[[1]]
    pat_tokens <- pat_tokens[nzchar(pat_tokens)]
    if (!length(pat_tokens)) next
    n <- length(tokens); m <- length(pat_tokens); if (n < m) next
    starts_idx <- which(sapply(seq_len(n - m + 1L),
                               function(s) all(tokens[s:(s + m - 1L)] == pat_tokens)))
    if (!length(starts_idx)) next
    ends_idx <- starts_idx + m - 1L
    as_bed <- rbind(as_bed,
                    data.frame(V1 = pat_str, V2 = starts_idx, V3 = ends_idx, stringsAsFactors = FALSE))
  }
  if (!nrow(as_bed)) return(as_bed)
  as_bed$pattern_start <- starts[as_bed$V2]
  as_bed$pattern_end   <- ends[as_bed$V3]
  as_bed
}

# ---------- process all CSVs ----------
data_files <- dir(CSV_DIR, pattern = "^chr[0-9]+_[0-9]+_bins_nam_groups\\.csv$", full.names = TRUE)
if (!length(data_files)) {
  message("No *_bins_nam_groups.csv found under ", CSV_DIR)
  quit(save = "no", status = 0)
}

for (csvf in data_files) {
  dat <- read.csv(csvf, stringsAsFactors = FALSE)
  if (!all(c("nam", "clust") %in% names(dat))) {
    warning("Skipping (missing nam/clust): ", csvf); next
  }
  pattern_out <- pattern_string_v2(dat)

  HOR_patterns_all <- data.frame(
    chr=character(), array_start=character(), array_end=character(),
    bin_start=character(), bin_end=character(),
    kmer=character(), Freq=integer(), charuni=integer(),
    char_prop_m=numeric(), char_num=integer(),
    stringsAsFactors = FALSE
  )
  HOR_patterns_all_bed <- data.frame(
    chr=character(), array_start=character(), array_end=character(),
    bin_start=character(), bin_end=character(),
    HOR_pattern=character(), start_mon=integer(), end_mon=integer(),
    pattern_start=numeric(), pattern_end=numeric(),
    stringsAsFactors = FALSE
  )
  string_out <- data.frame(chr=character(), array_start=character(), array_end=character(),
                           bin_start=character(), bin_end=character(),
                           threshold=character(), pattern=character(),
                           stringsAsFactors = FALSE)

  for (p in unique(pattern_out$id)) {
    sub <- pattern_out[pattern_out$id %in% p, ]

    letter_path <- file.path(LETTER_DIR, paste0(p, ".out"))
    orig_dat <- read_letters_2col(letter_path)
    if (is.null(orig_dat)) { warning("Missing/unparseable letter: ", letter_path); next }

    # coerce join keys
    sub <- sub %>% mutate(chr=as.character(chr),
                          array_start=as.character(array_start),
                          array_end=as.character(array_end),
                          letter=as.character(letter))
    orig_dat <- orig_dat %>% mutate(chr=as.character(chr),
                                    array_start=as.character(array_start),
                                    array_end=as.character(array_end),
                                    letter=as.character(letter),
                                    monomer_start=as.numeric(monomer_start),
                                    monomer_end=as.numeric(monomer_end))

    orig_dat_M <- merge(orig_dat, sub,
                        by = c("chr","array_start","array_end","letter"),
                        all = TRUE)
    orig_dat_M$new_letter[is.na(orig_dat_M$new_letter)] <- "Z"
    final_bin <- orig_dat_M %>% arrange(monomer_start)

    pattern <- paste(final_bin$new_letter, collapse = " ")

    new <- data.frame(chr=NA, array_start=NA, array_end=NA,
                      bin_start=NA, bin_end=NA, threshold=NA, pattern=NA,
                      stringsAsFactors = FALSE)
    p_sep <- stringr::str_split_fixed(p, "_", 5)
    new$chr         <- p_sep[,1]
    new$array_start <- p_sep[,2]
    new$array_end   <- p_sep[,3]
    new$bin_start   <- p_sep[,4]
    new$bin_end     <- p_sep[,5]
    new$threshold   <- suppressWarnings(max(as.numeric(sub$thr), na.rm = TRUE))
    new$pattern     <- pattern
    string_out <- rbind(string_out, new)

    HOR_candidates <- filtered_kmer_counts(pattern)
    HOR_patterns   <- filtered_kmer_counts_2(HOR_candidates)

    if (!is.null(HOR_patterns) && nrow(HOR_patterns) > 0) {
      HOR_patterns_all <- rbind(
        HOR_patterns_all,
        data.frame(
          chr = new$chr, array_start = new$array_start, array_end = new$array_end,
          bin_start = new$bin_start, bin_end = new$bin_end,
          kmer = HOR_patterns$out,
          Freq = HOR_patterns$Freq,
          charuni = HOR_patterns$charuni,
          char_prop_m = HOR_patterns$char_prop_m,
          char_num = HOR_patterns$char_num,
          stringsAsFactors = FALSE
        )
      )

      bed_df <- bed_convert(
        data.frame(HOR_patterns$out, stringsAsFactors = FALSE),
        pattern,
        starts = final_bin$monomer_start,
        ends   = final_bin$monomer_end
      )
      if (nrow(bed_df) > 0) {
        HOR_patterns_all_bed <- rbind(
          HOR_patterns_all_bed,
          data.frame(
            chr = new$chr, array_start = new$array_start, array_end = new$array_end,
            bin_start = new$bin_start, bin_end = new$bin_end,
            HOR_pattern = bed_df$V1,
            start_mon   = bed_df$V2, end_mon = bed_df$V3,
            pattern_start = bed_df$pattern_start,
            pattern_end   = bed_df$pattern_end,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  base <- sub("\\.csv$", "", basename(csvf))
  write.csv(string_out,       file.path(STRING_DIR, paste0(base, "_patterns_per_bin.csv")), row.names = FALSE)
  write.csv(HOR_patterns_all, file.path(STRING_DIR, paste0(base, "_HOR_patterns.csv")),       row.names = FALSE)
  write.csv(HOR_patterns_all_bed,
            file.path(STRING_DIR, paste0(base, "_HOR_patterns_bed.csv")),                     row.names = FALSE)

  message("[✓] Wrote: ", base)
}