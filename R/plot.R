#' Visualize a matrix as a heatmap and save it to a file
#'
#' Creates a heatmap visualization of a numeric matrix with a jet-like color palette.
#'
#' @param z A numeric matrix to be plotted. Must be numeric and two-dimensional.
#' @param filepath Character. The file path where the plot image will be saved if `save.image` is `TRUE`. Default is "plot.tiff" in the working directory.
#' @param save.image Logical. If `TRUE`, saves the plot to a file specified by `filepath`. Default is `FALSE`.
#' @param width Integer. Width of the saved image, in pixels (for raster formats) or inches (for vector formats). Default is 2500.
#' @param height Integer. Height of the saved image, in pixels (for raster formats) or inches (for vector formats). Default is 2350.
#' @param res Integer. Resolution of the saved image, in ppi (for raster formats only). Default is 300.
#' @param format Character. Format of the saved image. Options are `"tiff"`, `"png"`, `"jpeg"`, `"pdf"`, `"eps"`, or `"svg"`. Default is `"tiff"`.
#' @param colorbar.range Numeric vector of length 2. Specifies the minimum and maximum values for the color bar range. If `NULL` (default), the range is automatically determined from the input matrix `z`.
#' @param ... Additional graphical parameters to pass to the `image` function, such as axis properties.
#'
#' @details
#' This function plots a heatmap of a matrix with a color palette similar to the jet color map in MATLAB.
#' If `save.image` is `TRUE`, the plot is saved in the format and location specified by `filepath` and `format`.
#' The `colorbar.range` parameter allows users to manually set the range of the color bar, which is useful for ensuring consistent scaling across multiple plots.
#'
#' @return The function does not return any values. If `save.image` is `TRUE`, a message is printed to the console indicating the file path where the image was saved.
#'
#' @importFrom grDevices tiff png jpeg svg colorRampPalette dev.off pdf postscript
#' @importFrom graphics axis box layout par plot segments image
#' @export
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10)
#' plotMatrix(mat)

plotMatrix <- function(z,
                       filepath = file.path(getwd(), "plot.tiff"),
                       save.image = FALSE,
                       width = 2500,
                       height = 2350,
                       res = 300,
                       format = "tiff",
                       colorbar.range = NULL,  # New parameter for color bar range
                       ...) {
  # Sanity check
  if (!is.matrix(z)) {
    stop("Input must be a matrix.")
  }

  if (!is.numeric(z)) {
    stop("Data must be numeric.")
  }

  # Define values for plotting
  x <- seq(ncol(z))
  y <- seq(nrow(z))

  # Set zlim based on colorbar.range or the range of z
  if (!is.null(colorbar.range)) {
    if (!is.numeric(colorbar.range) || length(colorbar.range) != 2) {
      stop("colorbar.range must be a numeric vector of length 2 (c(min, max)).")
    }
    zlim <- colorbar.range
  } else {
    zlim <- range(z)
  }

  if (save.image) {
    # Start the appropriate graphics device
    switch(tolower(format),
           "tiff" = tiff(filename = filepath, width = width, height = height, res = res, compression = "lzw"),
           "png" = png(filename = filepath, width = width, height = height, res = res),
           "jpeg" = jpeg(filename = filepath, width = width, height = height, res = res),
           "svg" = svg(filename = filepath, width = width / 250, height = height / 250),
           "eps" = postscript(file = filepath, width = width / 250, height = height / 250, paper = "special", horizontal = FALSE),
           "pdf" = pdf(file = filepath, width = width / 250, height = height / 250),
           stop("Unsupported file format. Use 'tiff', 'png', 'jpeg', 'svg', 'eps', or 'pdf'.")
    )
  }

  # Create a jet-like palette
  color_palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(256)
  color_legend <- data.frame("color" = color_palette,
                             "value" = seq(from = zlim[1], to = zlim[2], along.with = 1:256))

  # Main plot
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.5, 1))
  par(mar = c(5, 4, 4, 1), ...)

  image(x = x,
        y = y,
        z = t(z[nrow(z):1, ]),
        zlim = zlim,
        axes = FALSE,
        ann = FALSE,
        col = color_palette,
        ...)

  # Handle xy-axis ticks
  axp.x <- if ("xaxp" %in% names(list(...))) {
    list(...)[["xaxp"]]  # Use custom x-axis parameters if provided
  } else {
    par("xaxp")  # Otherwise, use default parameters
  }
  at.x <- pretty(axp.x[1]:axp.x[2], axp.x[3])  # Generate pretty tick marks
  axis(1, at.x)  # Draw x-axis with pretty ticks

  axp.y <- if ("yaxp" %in% names(list(...))) {
    list(...)[["yaxp"]]
  } else {
    par("yaxp")
  }
  at.y <- pretty(axp.y[1]:axp.y[2], axp.y[3])
  at.y.r <- max(ncol(z)) - at.y + 1
  axis(2, at.y.r, labels = at.y, las = 1)

  box()

  # Add the legend
  par(mar = c(5, 0, 4, 5))
  plot(x = rep(1, length(color_legend$value)),
       y = color_legend$value,
       xlim = c(0, 1),
       col = color_legend$color,
       type = "n", xaxs = "i", yaxs = "i",
       ann = FALSE, axes = FALSE)

  segments(x0 = 0, x1 = 1,
           y0 = color_legend$value, y1 = color_legend$value,
           col = color_legend$color, lwd = 5)

  at.col.min <- ceiling(zlim[1] * 10) / 10  # Round up to nearest 0.1
  at.col.max <- floor(zlim[2] * 10) / 10      # Round down to nearest 0.1
  at.col <- seq(at.col.min, at.col.max, length.out = 5)

  axis(side = 4, at.col, lwd = 0, las = 2, line = -0.75)

  box()

  invisible(NULL)

  if (save.image) {
    dev.off()
    cat("Figure has been saved to:", normalizePath(filepath), "\n")
  }
}
