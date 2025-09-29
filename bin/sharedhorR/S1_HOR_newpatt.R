#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(stringr)
  library(igraph)
  library(caret)
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
IN_DIR  <- argval("--input",  default = ".")
OUT_DIR <- argval("--outdir", default = NULL)
if (is.null(OUT_DIR) || OUT_DIR == "") stop("Missing --outdir")
IN_DIR  <- normalizePath(IN_DIR,  mustWork = TRUE)
OUT_DIR <- normalizePath(OUT_DIR, mustWork = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- I/O dirs ----------
plot_dir <- file.path(OUT_DIR, "network_share_HOR")
csv_dir  <- file.path(OUT_DIR, "csv_share_HOR")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir,  showWarnings = FALSE, recursive = TRUE)

# ---------- Inputs: BLAT sub summaries ----------
data_files <- dir(IN_DIR,
                  pattern = "chr[0-9]+\\.cons\\.(9[0-9])\\.blat\\.sub$",
                  recursive = TRUE, full.names = TRUE)

if (!length(data_files)) {
  message("No .blat.sub files found in: ", IN_DIR)
  quit(save = "no", status = 0)
}

for (file in data_files) {
  bx  <- basename(file)
  thr <- as.numeric(sub(".*\\.cons\\.([0-9]+)\\.blat.*", "\\1", bx)) / 100

  # read pairwise data
  df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
  colnames(df) <- c("match", "nam1", "len1", "nam2", "len2")

  # jaccard
  df <- df %>% mutate(
    match = as.numeric(match),
    len1  = as.numeric(len1),
    len2  = as.numeric(len2),
    jacc  = match / (len1 + len2 - match)
  )

  # max jacc per pair (undirected)
  df2      <- aggregate(jacc ~ nam1 + nam2, df, max)
  rev_df2  <- df2; colnames(rev_df2) <- c("nam2", "nam1", "jacc")
  df2_full <- rbind(df2, rev_df2)

  # jaccard matrix
  jac_mat <- acast(df2_full, nam1 ~ nam2, value.var = "jacc", fun.aggregate = max, fill = 0)
  diag(jac_mat) <- 1

  # de-duplicate highly correlated
  dup  <- character(0)
  keep <- colnames(jac_mat)
  if (nrow(jac_mat) > 1) {
    dup  <- tryCatch(findCorrelation(jac_mat, cutoff = 0.9999999, names = TRUE),
                     error = function(e) character(0))
    keep <- setdiff(colnames(jac_mat), dup)
  }
  nams <- data.frame(nam = keep, stringsAsFactors = FALSE)
  nams$count <- NA_real_

  # count identicals (val==1) in original matrix
  if (length(keep) > 1) {
    for (j in seq_len(nrow(nams))) {
      row <- as.data.frame(t(jac_mat[nams$nam[j], , drop = FALSE]))
      colnames(row) <- "val"
      nams$count[j] <- sum(row$val == 1, na.rm = TRUE)
    }
  } else if (length(keep) == 1) {
    nams$count[1] <- sum(jac_mat[keep, ] == 1, na.rm = TRUE)
  }

  # thresholded adjacency
  jac_mat_sub <- jac_mat[keep, keep, drop = FALSE]
  jac_mat_sub[jac_mat_sub <  thr] <- 0
  jac_mat_sub[jac_mat_sub >= thr] <- 1

  network <- graph_from_adjacency_matrix(jac_mat_sub, mode = "undirected",
                                         weighted = TRUE, diag = FALSE)

  # components
  comp <- components(network)
  group_df <- data.frame(nam = names(comp$membership),
                         cluster = comp$membership,
                         stringsAsFactors = FALSE)

  # cluster counts (bring identical-count if you need)
  group_counts <- merge(group_df, nams, by = "nam", all.x = TRUE)

  # lookup for duplicates → inherit cluster of best neighbor
  cluster_lookup <- setNames(rep(NA_integer_, length(colnames(jac_mat))), colnames(jac_mat))
  if (nrow(group_df) > 0) cluster_lookup[group_df$nam] <- group_df$cluster
  if (length(dup) > 0) {
    for (d in dup) {
      simvals <- jac_mat[d, keep]
      best <- keep[which.max(simvals)]
      cluster_lookup[d] <- cluster_lookup[best]
    }
  }

  # expand seed groups including identicals
  nams_group <- data.frame(val = integer(0),
                           nam = character(0),
                           group = character(0),
                           clust = integer(0),
                           stringsAsFactors = FALSE)

  if (length(keep) > 1) {
    for (k in seq_len(nrow(nams))) {
      seed <- nams$nam[k]
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

    nfile <- gsub("^([^\\.]+).*\\.(\\d+)\\.blat\\.sub$", "\\1_\\2", bx)
    out_csv <- file.path(csv_dir, paste0(nfile, "_bins_nam_groups.csv"))
    write.csv(nams_group, out_csv, row.names = FALSE, quote = FALSE)
  }

  # plotting
  lay <- layout_with_fr(network)
  lay_rescaled <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
  idx <- match(V(network)$name, nams$nam)
  vcol <- ifelse(!is.na(nams$count[idx]) & nams$count[idx] > 1, "red", "gray70")
  V(network)$color <- vcol
  V(network)$frame.color <- "black"
  V(network)$label <- 1:vcount(network)
  V(network)$label.color <- "black"
  V(network)$label.cex <- 0.8

  pdf(file.path(plot_dir, paste0(nfile, "_network_plot.pdf")), width = 10, height = 10)
  plot(network,
       layout = lay_rescaled,
       rescale = FALSE,
       vertex.size = 8,
       vertex.label = V(network)$label,
       vertex.label.cex = 0.7,
       vertex.label.color = "black")
  mtext(paste("Threshold", thr), side = 3, line = 2, adj = 1, col = "black", cex = 1)
  dev.off()
}