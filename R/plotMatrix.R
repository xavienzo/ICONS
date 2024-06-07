#' Visualize a matrix as a heatmap and save it to a file
#'
#' This function creates a heatmap to visualize a numeric matrix. The heatmap is saved to a file
#' at the specified filepath in one of several supported formats.
#'
#' @param data A numeric matrix.
#' @param filepath The output file path. 
#' @param width Width of the output image in pixels, default is 850.
#' @param height Height of the output image in pixels, default is 800.
#' @param palette Type of color palette to use: "jet" or "viridis", default is "viridis".
#' @param format File format for saving the image, default is "tiff". Supported formats: "tiff", "png", "jpeg", "svg".
#' @return Nothing is explicitly returned; an image file is saved to the designated filepath.
#' @importFrom grDevices colorRampPalette dev.off jpeg png svg tiff
#' @importFrom graphics axis image mtext par segments
#' @importFrom stats cov quantile
#' @export
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10)
#' plotMatrix(mat, "heatmap.tiff")
plotMatrix <- function(data,
                       filepath,
                       width = 850,
                       height = 800,
                       palette = "viridis",
                       format = "tiff"
                       ) {
  # Sanity check
  if (!is.matrix(data)) {
    stop("Input must be a matrix.")
  }
  
  if (!is.numeric(data)) {
    stop("Data must be numeric.")
  }
  
  # Define color palettes
  if (palette == "jet") {
    col_pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"))(256)
  } else if (palette == "viridis") {
    col_pal <- colorRampPalette(c("#440154", "#21908C", "#FDE725"))(256)
  } else {
    stop("Unknown palette type. Use 'jet' or 'viridis'.")
  }
  # Record color to create legend
  col.key <- data.frame("color" = col_pal,
                        "value" = seq(from = min(data), to = max(data), along.with = 1:256))
  # Start the appropriate graphics device
  switch(tolower(format),
         "tiff" = tiff(filename = filepath, width = width, height = height, compression = "lzw"),
         "png" = png(filename = filepath, width = width, height = height, compression = "lzw"),
         "jpeg" = jpeg(filename = filepath, width = width, height = height),
         #"pdf" = pdf(filename = filepath, width = width, height = height),
         "svg" = svg(filename = filepath, width = width, height = height),
         stop("Unsupported file format. Use 'tiff', 'png', 'jpeg', or 'svg'.")
  )
  
  # Create the figure, leaving margins for legend
  par(fig = c(0,.9,0,1), mar = c(2,2,2,0))
  image(data, col = col_pal, useRaster = TRUE, axes = F, ann = F)
  # Add the legend
  par(fig = c(.9,1,.3,.7), mar = c(1,1,1,2.5), new = T)
  plot(x = rep(1,length(col.key$value)), y = col.key$value, xlim = c(0,1), col = col.key$color, type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
  segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
  axis(side = 4,lwd = 0, las = 2, line = -.75)
  mtext(text = "", adj = 0, line = 1, cex = 1.2)
  dev.off()
}