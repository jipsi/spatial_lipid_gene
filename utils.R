# utils.R
# Shared utility functions for the spatial_lipid_gene analysis pipeline.

#' Find the valley between the two highest density peaks (positive values only).
#'
#' @param data Numeric vector.
#' @return A list with valley_x, density, peak1_x, peak2_x, or NULL.
find_main_valley <- function(data) {
  data <- data[data > 0]

  if (length(data) < 10) {
    warning("Not enough positive values in dataset")
    return(NULL)
  }

  dens <- density(data, from = 0)

  derivative <- diff(dens$y)
  peaks <- which(diff(sign(derivative)) == -2) + 1

  if (length(peaks) >= 2) {
    peak_heights <- dens$y[peaks]
    top_two_peaks <- peaks[order(peak_heights, decreasing = TRUE)[1:2]]
    top_two_peaks <- sort(top_two_peaks)

    valleys <- which(diff(sign(derivative)) == 2) + 1
    valley_between <- valleys[valleys > top_two_peaks[1] &
                                valleys < top_two_peaks[2]]

    if (length(valley_between) > 0) {
      valley_depths <- dens$y[valley_between]
      main_valley <- valley_between[which.min(valley_depths)]

      return(list(
        valley_x = dens$x[main_valley],
        density = dens,
        peak1_x = dens$x[top_two_peaks[1]],
        peak2_x = dens$x[top_two_peaks[2]]
      ))
    }
  }
  return(NULL)
}


#' Plot a histogram with density curve, peaks, and the valley between them.
#'
#' @param data Numeric vector (positive values only).
#' @param title Character string for the plot title.
#' @return The x-coordinate of the valley, or NULL.
plot_with_valley <- function(data, title) {
  data <- data[data > 0]
  if (length(data) < 10) {
    stop("Not enough positive values in dataset")
  }

  result <- find_main_valley(data)

  if (!is.null(result)) {
    hist(data,
         prob = TRUE,
         breaks = "FD",
         col = "lightblue",
         border = "white",
         main = title,
         xlab = "Value (>0)",
         ylab = "Density",
         xlim = c(0, max(data)))

    lines(result$density, col = "darkblue", lwd = 2)
    abline(v = result$valley_x, col = "red", lwd = 2, lty = 2)
    abline(v = result$peak1_x, col = "darkgreen", lwd = 2, lty = 2)
    abline(v = result$peak2_x, col = "darkgreen", lwd = 2, lty = 2)

    text(result$valley_x,
         max(result$density$y) / 2,
         sprintf("Valley\n%.2f", result$valley_x),
         pos = 4, col = "red")

    grid(lty = "dotted")

    legend("topright",
           legend = c("Density", "Valley", "Peaks"),
           col = c("darkblue", "red", "darkgreen"),
           lwd = 2,
           lty = c(1, 2, 2))

    cat("Summary of positive values:\n")
    print(summary(data))
    cat("\nNumber of positive values:", length(data), "\n")

    return(result$valley_x)
  } else {
    warning("Could not find clear valley between two peaks")
    return(NULL)
  }
}


#' Compute Spearman correlations between a gene and all other genes.
#'
#' @param obj A Seurat object with an SCT assay.
#' @param gene_name Character; the gene to correlate against all others.
#' @return A list with two data.frames: corr_<gene_name> and acorr_<gene_name>.
fn_get_corr_mat <- function(obj, gene_name) {
  mat_count <- as.matrix(obj[['SCT']]@data)
  type <- "spearman"
  count_gene <- mat_count[gene_name, ]

  correlation_mat <- matrix(nrow = nrow(mat_count), ncol = 2)
  rownames(correlation_mat) <- rownames(mat_count)

  for (row in seq_len(nrow(mat_count))) {
    correlation <- stats::cor.test(count_gene, mat_count[row, ], method = type)
    correlation_mat[row, 1] <- correlation$estimate
    correlation_mat[row, 2] <- correlation$p.value
  }

  correlation_mat <- as.data.frame(correlation_mat)
  colnames(correlation_mat) <- c("corr_estimate", "pvalue")
  correlation_mat$gene <- rownames(correlation_mat)

  correlation_df <- sqldf::sqldf("SELECT gene, corr_estimate, pvalue
                                  FROM correlation_mat
                                  WHERE pvalue < 0.05
                                  AND corr_estimate > 0.1
                                  ORDER BY corr_estimate DESC")

  anti_correlation_df <- sqldf::sqldf("SELECT gene, corr_estimate, pvalue
                                       FROM correlation_mat
                                       WHERE pvalue < 0.05
                                       AND corr_estimate < -0.1
                                       ORDER BY corr_estimate ASC")

  list_corr_mat <- list()
  list_corr_mat[[paste0("corr_", gene_name)]] <- correlation_df
  list_corr_mat[[paste0("acorr_", gene_name)]] <- anti_correlation_df
  return(list_corr_mat)
}
