#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)   # for acast
  library(stringr)
  library(igraph)
  library(caret)      # for findCorrelation
})

`%!in%` <- Negate(`%in%`)
options(scipen = 999)

# ---------------- CLI ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) return(TRUE)
  args[i + 1]
}

INPUT_DIR   <- get_arg("--input")
BINS_FILE   <- get_arg("--bins")
OUTDIR      <- get_arg("--outdir")
DEFAULT_THR <- as.numeric(get_arg("--default-thr", "0.35"))

if (is.null(INPUT_DIR) || is.null(BINS_FILE) || is.null(OUTDIR)) {
  stop(paste(
    "Usage: S3_network_summary_HOR_mergedbin.R",
    "--input DIR --bins FILE --outdir DIR [--default-thr 0.35]"
  ), call. = FALSE)
}

if (!dir.exists(INPUT_DIR)) stop("Input dir not found: ", INPUT_DIR)
if (!file.exists(BINS_FILE)) stop("Bins file not found: ", BINS_FILE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

plot_dir <- file.path(OUTDIR, "network_plots_mergebin")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("[*] INPUT  : ", normalizePath(INPUT_DIR))
message("[*] BINS   : ", normalizePath(BINS_FILE))
message("[*] OUTDIR : ", normalizePath(OUTDIR))
message("[*] Default threshold if missing in bins: ", DEFAULT_THR)

# ---------------- Inputs ----------------
data_files <- dir(INPUT_DIR, pattern = "\\.blat\\.sub$", recursive = TRUE, full.names = TRUE)
if (!length(data_files)) stop("No *.blat.sub found under: ", INPUT_DIR)

bins_fil <- read.table(BINS_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)

# Expect these columns in bins file
need_bins_cols <- c("chr", "start", "end", "array", "threshold")
miss_bins <- setdiff(need_bins_cols, colnames(bins_fil))
if (length(miss_bins)) {
  stop("Bins file missing required columns: ", paste(miss_bins, collapse = ", "))
}

info_summary_all  <- list()
label_summary_all <- list()

# ---------------- Main loop ----------------
for (i in seq_along(data_files)) {
  file <- data_files[i]
  bx   <- basename(file)

  # Parse metadata (works for "chr..." and "all_chr...")
  m <- stringr::str_match(
    bx,
    "^(?:all_)?(chr[^_]+)_(\\d+)_(\\d+)_(\\d+)_(\\d+)\\.blat\\.sub$"
  )
  if (any(is.na(m))) stop("Filename does not match expected pattern: ", bx)

  chr         <- m[2]
  array_start <- as.numeric(m[3])
  array_end   <- as.numeric(m[4])
  bin_start   <- as.numeric(m[5])
  bin_end     <- as.numeric(m[6])
  array_str   <- paste0(array_start, "-", array_end)

  # Clean stem for plots/labels
  sfile_base <- sub("\\.blat\\.sub$", "", bx)  # drop suffix
  sfile      <- sub("^all_", "", sfile_base)   # drop optional 'all_'

  # Threshold from bins (fallback to DEFAULT_THR)
  df_thr <- data.frame(array = array_str, chr = chr, start = bin_start, end = bin_end,
                       stringsAsFactors = FALSE)
  df_joined <- df_thr %>%
    left_join(bins_fil, by = c("chr"="chr","start"="start","end"="end","array"="array")) %>%
    select(array, chr, start, end, threshold)
  thr <- df_joined$threshold[1]
  if (is.na(thr) || !is.finite(thr)) thr <- DEFAULT_THR

  # Read BLAT-sub (match, nam1, len1, nam2, len2)
  df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
  if (ncol(df) < 5) stop("Bad .blat.sub format (need 5 cols) in: ", bx)
  colnames(df) <- c("match", "nam1", "len1", "nam2", "len2")

  df$nam1 <- gsub("[:\\-]", "_", df$nam1)
  df$nam2 <- gsub("[:\\-]", "_", df$nam2)

  df <- df %>% mutate(
    match = as.numeric(match),
    len1  = as.numeric(len1),
    len2  = as.numeric(len2),
    jacc  = match / (len1 + len2 - match)
  )

  # Max jacc per pair (undirected)
  df2      <- aggregate(jacc ~ nam1 + nam2, df, max)
  rev_df2  <- df2; colnames(rev_df2) <- c("nam2","nam1","jacc")
  df2_full <- rbind(df2, rev_df2)

  # Jaccard matrix
  jac_mat <- acast(df2_full, nam1 ~ nam2, value.var = "jacc", fun.aggregate = max, fill = 0)
  diag(jac_mat) <- 1

  # Deduplicate near-identicals
  dup  <- character(0)
  keep <- colnames(jac_mat)
  if (nrow(jac_mat) > 1) {
    dup  <- tryCatch(findCorrelation(jac_mat, cutoff = 0.9999999, names = TRUE),
                     error = function(e) character(0))
    keep <- setdiff(colnames(jac_mat), dup)
  }
  nams <- data.frame(nam = keep, stringsAsFactors = FALSE)
  nams$count <- NA_real_

  if (length(keep) > 1) {
    for (j in seq_len(nrow(nams))) {
      row <- as.data.frame(t(jac_mat[nams$nam[j], , drop = FALSE]))
      colnames(row) <- "val"
      nams$count[j] <- sum(row$val == 1, na.rm = TRUE)
    }
  } else if (length(keep) == 1) {
    nams$count[1] <- sum(jac_mat[keep, ] == 1, na.rm = TRUE)
  }

  # ---- Trivial case: <=1 kept monomer ----
  if (length(keep) <= 1) {
    num_monomers <- length(unique(df$nam1))

    info <- data.frame(
      array        = array_str,
      chr          = chr,
      start        = bin_start,
      end          = bin_end,
      threshold    = thr,
      num_monomers = num_monomers,
      dup_monomers = 0,
      uncon_clust  = 1,
      prop_clust1  = 1,
      prop_clust2  = 0,
      modularity   = 0,
      mean_dist    = 0,
      avg_jacc     = if (nrow(df2_full) > 0) mean(df2_full$jacc) else 1,
      max_contract = 1,
      stringsAsFactors = FALSE
    )
    info_summary_all[[length(info_summary_all) + 1]] <- info

    label_df <- data.frame(
      bin       = sfile,   # clean stem
      label     = "A",
      threshold = thr,
      stringsAsFactors = FALSE
    )
    label_summary_all[[length(label_summary_all) + 1]] <- label_df

    next
  }

  # ---- Build graph at threshold thr ----
  jac_mat_sub <- jac_mat[keep, keep, drop = FALSE]
  jac_mat_sub[jac_mat_sub <  thr] <- 0
  jac_mat_sub[jac_mat_sub >= thr] <- 1

  network <- graph_from_adjacency_matrix(jac_mat_sub, mode = "undirected", weighted = TRUE, diag = FALSE)

  comp <- components(network)
  group_df <- data.frame(nam = names(comp$membership), cluster = comp$membership, stringsAsFactors = FALSE)

  # rank clusters by size
  cluster_order <- group_df %>%
    count(cluster, name = "size") %>%
    arrange(desc(size)) %>%
    mutate(new_cluster = row_number())

  group_df <- group_df %>%
    left_join(cluster_order[, c("cluster","size","new_cluster")], by = "cluster") %>%
    mutate(cluster = new_cluster) %>%
    select(nam, cluster)

  kmax <- if (nrow(group_df) > 0) max(group_df$cluster) else 0
  cluster_labels <- c(LETTERS[1:25], 1:9, letters[1:25])
  cluster_labels <- if (kmax <= length(cluster_labels)) cluster_labels[1:kmax] else paste0("C", seq_len(kmax))

  group_counts <- merge(group_df, nams, by = "nam", all.x = TRUE)
  cluster_sum <- group_counts %>%
    group_by(cluster) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    arrange(cluster)

  if (kmax > 0) {
    cluster_sum$label <- cluster_labels[cluster_sum$cluster]
    group_counts <- merge(group_counts, cluster_sum[, c("cluster","label")], by = "cluster", all.x = TRUE)
  } else {
    group_counts$label <- character(0)
  }

  # duplicates inherit best kept neighbor
  label_map <- group_counts[, c("nam","label")]
  if (length(dup) > 0) {
    dup_df <- data.frame(nam = dup, stringsAsFactors = FALSE)
    dup_df$label <- sapply(dup, function(d) {
      simvals <- jac_mat[d, keep]
      best <- keep[which.max(simvals)]
      label_map$label[label_map$nam == best]
    })
    group_counts_full <- rbind(group_counts[, c("nam","label")], dup_df)
  } else {
    group_counts_full <- group_counts[, c("nam","label")]
  }

  # === Optional CSV like *_bins_nam_groups.csv ===
  cluster_lookup <- setNames(rep(NA_integer_, length(colnames(jac_mat))), colnames(jac_mat))
  if (nrow(group_df) > 0) cluster_lookup[group_df$nam] <- group_df$cluster
  if (length(dup) > 0) {
    for (d in dup) {
      simvals <- jac_mat[d, keep]
      best <- keep[which.max(simvals)]
      cluster_lookup[d] <- cluster_lookup[best]
    }
  }

  nams_group <- data.frame(val = integer(0), nam = character(0), group = character(0),
                           clust = integer(0), stringsAsFactors = FALSE)

  if (length(keep) > 1) {
    for (k in seq_len(nrow(group_counts))) {
      seed <- group_counts$nam[k]
      row_vals <- jac_mat[seed, , drop = FALSE]
      sub_nams <- names(which(row_vals[1, ] == 1))
      if (length(sub_nams) > 0) {
        sub_df <- data.frame(
          val   = 1L,
          nam   = sub_nams,
          group = seed,
          clust = as.integer(cluster_lookup[sub_nams]),
          stringsAsFactors = FALSE
        )
        nams_group <- rbind(nams_group, sub_df)
      }
    }
    nams_group$origin <- bx
    out_csv_dir <- file.path(OUTDIR, "out_csv")
    dir.create(out_csv_dir, showWarnings = FALSE)
    nfile   <- sub("^all_", "", sub("\\.blat\\.sub$", "", bx))
    out_csv <- file.path(out_csv_dir, paste0(nfile, "_bins_nam_groups.csv"))
    write.csv(nams_group, out_csv, row.names = FALSE, quote = FALSE)
  }

  # ---- plotting ----
  lay <- layout_with_fr(network)
  lay_rescaled <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1)

  idx <- match(V(network)$name, group_counts$nam)
  cnt <- group_counts$count[idx]
  vcol <- ifelse(!is.na(cnt) & cnt > 1, "red", "gray70")
  V(network)$color <- vcol
  V(network)$frame.color <- "black"
  V(network)$label <- 1:vcount(network)
  V(network)$label.color <- "black"
  V(network)$label.cex <- 0.8

  pdf(file.path(plot_dir, paste0(sfile, "thr", round(thr * 100), "_network_plot.pdf")),
      width = 10, height = 10)
  plot(network,
       layout = lay_rescaled,
       rescale = FALSE,
       vertex.size = 8,
       vertex.label = V(network)$label,
       vertex.label.cex = 0.7,
       vertex.label.color = "black")

  if (kmax > 0) {
    for (cl in seq_len(kmax)) {
      verts <- which(group_df$cluster == cl)
      if (length(verts) == 0) next
      vnames <- group_df$nam[verts]
      vind   <- match(vnames, V(network)$name)
      x <- max(lay_rescaled[vind, 1]); y <- max(lay_rescaled[vind, 2])
      text(x + 0.05, y + 0.05, labels = cluster_labels[cl], col = "blue", cex = 1.2, font = 2)
    }
  }

  group_counts_ord <- group_counts_full[order(group_counts_full$nam), ]
  label_string <- if (nrow(group_counts_ord) > 0) paste(group_counts_ord$label, collapse = "") else ""

  # Wrap the long label into multiple lines (tweak width/cex as you like)
  lab <- paste("Monomers", label_string)
  lab_split <- strwrap(lab, width = 60)  # adjust width
  for (ii in seq_along(lab_split)) {
    mtext(lab_split[ii], side = 3, line = ii, adj = 0, col = "black", cex = 0.5)
  }
  mtext(paste("Threshold", thr),
        side = 3, line = length(lab_split) + 0.5, adj = 1, col = "black", cex = 1)
  dev.off()

  # ---- summaries ----
  num_monomers <- length(unique(df$nam1))
  cluster_sizes <- if (nrow(cluster_sum) > 0) cluster_sum$count else numeric(0)
  cluster1 <- if (length(cluster_sizes) >= 1) cluster_sizes[1] else NA_real_
  cluster2 <- if (length(cluster_sizes) >= 2) cluster_sizes[2] else NA_real_

  dup_monomers <- (length(dup) + nrow(nams[!is.na(nams$count) & nams$count > 1, , drop = FALSE])) / num_monomers
  uncon_clust  <- if (num_monomers > 0) (kmax / num_monomers) else NA_real_
  prop_clust1  <- if (!is.na(cluster1)) cluster1 / num_monomers else 0
  prop_clust2  <- if (!is.na(cluster2)) cluster2 / num_monomers else 0
  modularity_v <- if (ecount(network) == 0 || kmax == 0) 0 else modularity(network, group_df$cluster)
  mean_dist_v  <- if (vcount(network) <= 1) 0 else mean_distance(network, directed = FALSE, unconnected = TRUE)
  avg_jacc_v   <- if (nrow(df2_full) > 0) mean(df2_full$jacc) else NA_real_
  max_contract <- if (num_monomers > 0) max(nams$count, na.rm = TRUE) / num_monomers else NA_real_

  info <- data.frame(
    array        = array_str,
    chr          = chr,
    start        = bin_start,
    end          = bin_end,
    threshold    = thr,
    num_monomers = num_monomers,
    dup_monomers = dup_monomers,
    uncon_clust  = uncon_clust,
    prop_clust1  = prop_clust1,
    prop_clust2  = prop_clust2,
    modularity   = modularity_v,
    mean_dist    = mean_dist_v,
    avg_jacc     = avg_jacc_v,
    max_contract = max_contract,
    stringsAsFactors = FALSE
  )
  info_summary_all[[length(info_summary_all) + 1]] <- info

  label_df <- data.frame(
    bin       = sfile,
    label     = label_string,
    threshold = thr,
    stringsAsFactors = FALSE
  )
  label_summary_all[[length(label_summary_all) + 1]] <- label_df
}

# ---------------- Write outputs ----------------
network_summary <- dplyr::bind_rows(info_summary_all)
label_summary   <- dplyr::bind_rows(label_summary_all)

out1 <- file.path(OUTDIR, "HOR_re-eval_clusters.txt")
out2 <- file.path(OUTDIR, "label_HOR_re-eval_clusters.txt")
write.table(network_summary, out1, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(label_summary,   out2, row.names = FALSE, sep = "\t", quote = FALSE)

message("[✓] Re-evaluation complete.")
message("    - ", normalizePath(out1))
message("    - ", normalizePath(out2))
message("    - Plots in: ", normalizePath(plot_dir))