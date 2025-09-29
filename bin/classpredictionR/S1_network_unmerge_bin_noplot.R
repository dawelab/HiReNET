#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(reshape2)
  library(caret)
  library(stringr)
  library(data.table)
})

options(scipen = 999)

# ---------------- CLI ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) return(TRUE) # boolean switch
  args[i + 1]
}

INPUT_DIR <- get_arg("--input")
OUTDIR    <- get_arg("--outdir")
PREFIX    <- get_arg("--prefix")
LDA_MODEL <- get_arg("--lda-model", default = NA_character_)
DO_PLOT   <- !is.null(get_arg("--plot", default = NULL))  # optional

if (is.null(INPUT_DIR) || is.null(OUTDIR) || is.null(PREFIX)) {
  stop("Usage: S1_network_unmerge_bin_noplot.R --input DIR --outdir DIR --prefix NAME [--lda-model model.rds] [--plot]\n",
       call. = FALSE)
}

if (!dir.exists(INPUT_DIR)) stop("Input dir not found: ", INPUT_DIR)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(OUTDIR, "network_plots")
if (DO_PLOT && !dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

message("[*] INPUT  : ", normalizePath(INPUT_DIR))
message("[*] OUTDIR : ", normalizePath(OUTDIR))
message("[*] PREFIX : ", PREFIX)
if (!is.na(LDA_MODEL)) {
  if (!file.exists(LDA_MODEL)) stop("LDA model not found: ", LDA_MODEL)
  message("[*] LDA    : ", normalizePath(LDA_MODEL))
} else {
  message("[*] LDA    : (none provided — prediction step will be skipped)")
}
message("[*] PLOT   : ", if (DO_PLOT) "yes" else "no")

# ------------- Inputs ---------------
data_files <- dir(INPUT_DIR, pattern = "\\.blat\\.sub$", recursive = TRUE, full.names = TRUE)
if (length(data_files) == 0) stop("No *.blat.sub files found under: ", INPUT_DIR)

# thresholds
clust_val <- c(.90, .91, .92, .93, .94, .95, .96, .97, .98, .99)

# per-threshold nested lists
info_summary_all <- vector("list", length(clust_val))  # each will be a list of data.frames

# ------------- Main loop ------------
for (file in data_files) {
  fname <- basename(file)
  sfile <- sub("\\.blat\\.sub$", "", fname)  # base for optional plots

    # skip zero-byte files early
  if (file.size(file) == 0) {
    message(sprintf("[skip] %s (empty file size)", basename(file)))
    next
  }

  # robust read
  df <- tryCatch(
    read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = ""),
    error = function(e) {
      message(sprintf("[skip] %s (read error: %s)", basename(file), e$message))
      return(NULL)
    }
  )
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("[skip] %s (no lines after read)", basename(file)))
    next
  }
  if (ncol(df) < 5) {
    message(sprintf("[skip] %s (expected 5 columns, found %d)", basename(file), ncol(df)))
    next
  }
  
  colnames(df) <- c("match", "nam1", "len1", "nam2", "len2")
  df$nam1 <- gsub("[:\\-]", "_", df$nam1)
  df$nam2 <- gsub("[:\\-]", "_", df$nam2)

  # Jaccard
  df <- df %>% mutate(jacc = match / (len1 + len2 - match))

  # max jacc per pair (undirected)
  df2      <- aggregate(jacc ~ nam1 + nam2, df, max)
  rev_df2  <- df2; colnames(rev_df2) <- c("nam2", "nam1", "jacc")
  df2_full <- rbind(df2, rev_df2)

  # matrix
  jac_mat <- acast(df2_full, nam1 ~ nam2, value.var = "jacc", fun.aggregate = max, fill = 0)
  diag(jac_mat) <- 1

  # deduplicate near-identicals
  dup  <- character(0)
  keep <- colnames(jac_mat)
  if (nrow(jac_mat) > 1) {
    dup  <- findCorrelation(jac_mat, cutoff = 0.9999999, names = TRUE)
    keep <- setdiff(colnames(jac_mat), dup)
    nams <- data.frame(nam = keep, stringsAsFactors = FALSE)
  } else {
    nams <- data.frame(nam = colnames(jac_mat), stringsAsFactors = FALSE)
  }
  nams$count <- NA_real_

  if (length(keep) > 1) {
    for (j in seq_len(nrow(nams))) {
      row <- as.data.frame(t(jac_mat[nams$nam[j], , drop = FALSE]))
      colnames(row) <- "val"
      row$nam <- rownames(row)
      sub <- row[row$val == 1, ]
      nams$count[j] <- nrow(sub)
    }
  }

  b <- basename(file)
  m <- regmatches(b, regexec("^(chr[^_]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)\\.blat\\.sub$", b))[[1]]
  if (length(m) == 0) {
    parts <- str_split(b, "_", simplify = TRUE)
    chr <- parts[1]
    array_start <- suppressWarnings(as.numeric(parts[2]))
    array_end   <- suppressWarnings(as.numeric(parts[3]))
    start <- suppressWarnings(as.numeric(parts[4]))
    end   <- suppressWarnings(as.numeric(gsub("\\.blat\\.sub$", "", parts[5])))
  } else {
    chr <- m[2]; array_start <- as.numeric(m[3]); array_end <- as.numeric(m[4]); start <- as.numeric(m[5]); end <- as.numeric(m[6])
  }
  array <- if (is.na(array_start) || is.na(array_end)) NA_character_ else paste0(array_start, "-", array_end)

  # trivial case (<=1 kept monomer)
  if (length(keep) <= 1) {
    num_monomers <- length(unique(df$nam1))
    for (k in seq_along(clust_val)) {
      thr <- clust_val[k]
      info <- data.frame(
        array        = array,
        chr          = chr,
        start        = start,
        end          = end,
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
      if (is.null(info_summary_all[[k]])) info_summary_all[[k]] <- list()
      info_summary_all[[k]][[length(info_summary_all[[k]]) + 1]] <- info
    }
    next
  }

  # normal case: iterate thresholds
  for (k in seq_along(clust_val)) {
    thr <- clust_val[k]

    jac_mat_sub <- jac_mat[keep, keep, drop = FALSE]
    jac_mat_sub[jac_mat_sub <  thr] <- 0
    jac_mat_sub[jac_mat_sub >= thr] <- 1

    network <- graph_from_adjacency_matrix(jac_mat_sub, mode = "undirected", weighted = TRUE, diag = FALSE)

    idx <- match(V(network)$name, nams$nam)
    vcol <- ifelse(nams$count[idx] > 1, "red", "gray70")
    V(network)$color <- vcol
    V(network)$frame.color <- "black"
    V(network)$label <- 1:vcount(network)
    V(network)$label.color <- "black"
    V(network)$label.cex <- 0.8

    comp <- components(network)
    group_df <- data.frame(nam = names(comp$membership), cluster = comp$membership, stringsAsFactors = FALSE)

    cluster_order <- group_df %>%
      count(cluster, name = "size") %>%
      arrange(desc(size)) %>%
      mutate(new_cluster = row_number())

    group_df <- group_df %>%
      left_join(cluster_order[, c("cluster", "size", "new_cluster")], by = "cluster") %>%
      mutate(cluster = new_cluster) %>%
      select(nam, cluster)

    kmax <- if (nrow(group_df) > 0) max(group_df$cluster) else 0

    group_counts <- merge(group_df, nams, by = "nam", all.x = TRUE)
    cluster_sum <- group_counts %>%
      group_by(cluster) %>%
      summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
      arrange(cluster)

    if (DO_PLOT) {
      lay <- layout_with_fr(network)
      lay_rescaled <- norm_coords(lay, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
      pdf(file.path(plot_dir, paste0(sfile, "_thr", round(thr * 100), "_network_plot.pdf")), width = 10, height = 10)
      plot(network,
           layout = lay_rescaled,
           rescale = FALSE,
           vertex.size = 8,
           vertex.label = V(network)$label,
           vertex.label.cex = 0.7,
           vertex.label.color = "black")
      mtext(paste("Threshold", thr), side = 3, line = 1.5, adj = 1, col = "black", cex = 1)
      dev.off()
    }

    # ---- Summaries ----
    num_monomers <- length(unique(df$nam1))
    cluster_sizes <- if (nrow(cluster_sum) > 0) cluster_sum$count else numeric(0)
    cluster1 <- if (length(cluster_sizes) >= 1) cluster_sizes[1] else NA_real_
    cluster2 <- if (length(cluster_sizes) >= 2) cluster_sizes[2] else NA_real_

    dup_monomers <- (length(dup) + nrow(nams[nams$count > 1, , drop = FALSE])) / num_monomers
    uncon_clust  <- if (num_monomers > 0) (kmax / num_monomers) else NA_real_
    prop_clust1  <- if (!is.na(cluster1)) cluster1 / num_monomers else 0
    prop_clust2  <- if (!is.na(cluster2)) cluster2 / num_monomers else 0

    modularity_v <- if (ecount(network) == 0 || kmax == 0) 0 else modularity(network, group_df$cluster)
    mean_dist_v  <- if (vcount(network) <= 1) 0 else mean_distance(network, directed = FALSE, unconnected = TRUE)
    avg_jacc_v   <- if (nrow(df2_full) > 0) mean(df2_full$jacc) else NA_real_
    max_contract <- max(nams$count, na.rm = TRUE) / num_monomers

    info <- data.frame(
      array        = array,
      chr          = chr,
      start        = start,
      end          = end,
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

    if (is.null(info_summary_all[[k]])) info_summary_all[[k]] <- list()
    info_summary_all[[k]][[length(info_summary_all[[k]]) + 1]] <- info
  }
}

# ---- Collapse and write summaries ----
network_summary_list <- lapply(info_summary_all, function(lst) if (length(lst)) dplyr::bind_rows(lst) else data.frame())

network_summary_combined <- dplyr::bind_rows(network_summary_list) %>%
  dplyr::arrange(chr, array, start, end, threshold)

write.table(network_summary_combined,
            file.path(OUTDIR, "network_summary_thresholds.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)

# ---- Optional: LDA prediction ----
if (!is.na(LDA_MODEL)) {
  fit2 <- readRDS(LDA_MODEL)

  network_summary_or_list <- list()
  for (i in seq_along(clust_val)) {
    network_summary <- network_summary_list[[i]]
    if (is.null(network_summary) || nrow(network_summary) == 0) next

    colnames(network_summary) <- c("array","chr","start","end","threshold","num_monomer",
                                   "dup_monomers","uncon_clust","prop_clust1","prop_clust2",
                                   "modularity","mean_dist","avg_jacc","max_contract")

    network_summary <- network_summary[network_summary$num_monomer >= 5, ]
    if (anyNA(network_summary$mean_dist))  network_summary$mean_dist[is.na(network_summary$mean_dist)] <- 10
    if (anyNA(network_summary$prop_clust2)) network_summary$prop_clust2[is.na(network_summary$prop_clust2)] <- 0
    if (anyNA(network_summary$modularity))  network_summary$modularity[is.na(network_summary$modularity)] <- 0

    pr <- predict(fit2, newdata = network_summary, type = "class")
    network_summary$lda_pred     <- pr$class
    network_summary$lda_pred_pos <- apply(pr$posterior, 1, max)

    network_summary_or_list[[i]] <- network_summary
  }

  if (length(network_summary_or_list)) {
    classes <- c("Exp","HOR","Order","Disorder")
    which_longest <- which.max(sapply(network_summary_or_list, nrow))
    nrows <- if (length(which_longest)) nrow(network_summary_or_list[[which_longest]]) else 0

    combined <- data.frame()
    for (i in seq_len(nrows)) {
      comp <- rbindlist(lapply(network_summary_or_list, function(x) if (nrow(x) >= i) x[i, ] else NULL))
      comp <- comp[comp$lda_pred_pos >= 0.80, ]
      if (nrow(comp) == 0) next
      idx <- match(comp$lda_pred, classes)
      comp$.prio <- ifelse(is.na(idx), length(classes) + 1, idx)
      comp <- comp[order(comp$.prio, -comp$lda_pred_pos), ]
      best <- comp[1, ]
      combined <- rbind(combined, best)
    }

    out_pred <- file.path(OUTDIR, paste0(PREFIX, "_bin_class.txt"))
    write.table(combined, out_pred, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

message("[✓] S1 complete. Outputs in: ", normalizePath(OUTDIR))