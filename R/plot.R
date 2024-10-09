#' Visualize a matrix as a heatmap and save it to a file
#'
#' This function creates a heatmap to visualize a numeric matrix. The heatmap is saved to a file
#' at the specified filepath in one of several supported formats.
#'
#' @param data A numeric matrix.
#' @param filepath The output file path, default is the current working directory
#' @param width Width of the output image in pixels, default is 850.
#' @param height Height of the output image in pixels, default is 800.
#' @param palette Type of color palette to use: "jet" or "viridis", default is "jet".
#' @param format File format for saving the image, default is "tiff". Supported formats: "tiff", "png", "jpeg", "svg".
#' @return Nothing is explicitly returned; an image file is saved to the designated file path.
#' @importFrom grDevices colorRampPalette dev.off jpeg png svg tiff
#' @importFrom graphics axis image mtext par segments
#' @importFrom stats cov quantile
#' @export
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10)
#' plotMatrix(mat)
plotMatrix <- function(data,
                       filepath = file.path(getwd(), "plot.tiff"),
                       width = 850,
                       height = 800,
                       palette = "jet",
                       format = "tiff"
) {
  # Sanity check ===================
  if (!is.matrix(data)) {
    stop("Input must be a matrix.")
  }

  if (!is.numeric(data)) {
    stop("Data must be numeric.")
  }

  # Define color palettes ===================
  if (palette == "jet") {
    color_palette <- colorRampPalette(c("blue", "cyan", "yellow", "red"))(256)
  } else if (palette == "viridis") {
    color_palette <- colorRampPalette(c("#440154", "#21908C", "#FDE725"))(256)
  } else {
    stop("Unknown palette type. Use 'jet' or 'viridis'.")
  }
  # Record color to create legend
  color_legend <- data.frame("color" = color_palette,
                             "value" = seq(from = min(data, na.rm = T), to = max(data, na.rm = T), along.with = 1:256))

  # Create the figure ===================

  # Start the appropriate graphics device
  switch(tolower(format),
         "tiff" = tiff(filename = filepath, width = width, height = height, compression = "lzw"),
         "png" = png(filename = filepath, width = width, height = height, compression = "lzw"),
         "jpeg" = jpeg(filename = filepath, width = width, height = height),
         "svg" = svg(filename = filepath, width = width, height = height),
         stop("Unsupported file format. Use 'tiff', 'png', 'jpeg', or 'svg'.")
  )

  # Leave margins for legend
  par(fig = c(0,.9,0,1), mar = c(2,2,2,0))
  image(t(data[nrow(data):1, ]), col = color_palette,
        useRaster = TRUE, axes = FALSE, ann = FALSE)

  # Function to select indices and normalize to [0, 1] range
  select_indices <- function(total_length, max_labels) {
    step <- max(1, ceiling(total_length / max_labels))
    indices <- seq(0, total_length, by = step)[-1]
    # Normalize indices to [0, 1]
    normalized_indices <- indices / total_length
    return(list(indices = indices, normalized = normalized_indices))
  }

  # Determine how many labels to show
  max_labels_x <- 5  # Maximum labels for x-axis
  max_labels_y <- 5  # Maximum labels for y-axis

  # Select indices
  x_labels <- select_indices(ncol(data), max_labels_x)
  y_labels <- select_indices(nrow(data), max_labels_y)

  # Add custom axis labels with normalized positions
  axis(1, at = x_labels$normalized, labels = x_labels$indices)  # X-axis
  axis(2, at = y_labels$normalized, labels = rev(y_labels$indices))  # Y-axis (reversed to match heatmap orientation)

  # Add the legend
  par(fig = c(.9,1,.2,.8), mar = c(1,1,1,2.5), new = TRUE)
  # par(fig = c(0.95, 1, 0, 1), mar = c(2, 0, 2, 2), new = TRUE)
  plot(x = rep(1,length(color_legend$value)),
       y = color_legend$value,
       xlim = c(0,1),
       col = color_legend$color,
       type = "n", xaxs = "i", yaxs = "i",
       ann = FALSE, axes = FALSE)
  segments(x0 = 0, x1 = 1,
           y0 = color_legend$value, y1 = color_legend$value,
           col = color_legend$color, lwd = 5)
  axis(side = 4, lwd = 0, las = 2, line = -.75)
  # mtext(text = "", adj = 0, line = 1, cex = 1.2)
  dev.off()
  cat("Figure has been saved to:", normalizePath(filepath), "\n")
}
