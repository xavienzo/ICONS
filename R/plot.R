#' Heat map of a matrix
#'
#' Draws a matrix the way `imagesc` does in MATLAB: row 1 at the top, a jet
#' colour ramp, and a colour bar on the right.  Optionally overlays the
#' community boundaries of a partition.
#'
#' Unlike ICONS 0.1.x this restores the graphics state on exit, so the device
#' is not left in a two-panel layout, and it downsamples very large matrices
#' rather than pushing millions of cells at a device that cannot resolve them.
#'
#' @param z A numeric matrix.  Need not be square.
#' @param partition Optional `"icons_partition"`.  When given, `z` is assumed
#'   to be already reordered by [reorder_matrix()] and community boundaries are
#'   drawn.
#' @param zlim Numeric length-2 colour-bar range.  Defaults to `range(z)`.
#'   Fix it across plots to make them comparable.
#' @param col Colour palette.  Defaults to a jet-like ramp of 256 colours.
#' @param max_cells Largest matrix side drawn at full resolution; larger
#'   matrices are block-averaged using the smallest whole-number block size that
#'   brings them under this. Set to `Inf` to disable.  Axis labels always refer
#'   to original variable indices, whether or not downsampling happened.
#' @param border Colour of the community boundary lines, or `NA` for none.
#' @param file Optional path to save to instead of drawing on the current
#'   device.  The extension picks the format unless `format` is given.
#' @param format One of `"tiff"`, `"png"`, `"jpeg"`, `"pdf"`, `"eps"`, `"svg"`.
#' @param width,height Size of the saved image, pixels for raster formats and
#'   inches for vector ones.
#' @param res Resolution in ppi for raster formats.
#' @param main,xlab,ylab Titles.
#' @param ... Passed to [graphics::image()].
#'
#' @return `NULL`, invisibly.  Called for the plot.
#'
#' @examples
#' data(sim)
#' W <- cor(sim)
#' part <- icons_detect(W, threshold = 0.6)
#'
#' plot_matrix(W, main = "Original")
#' plot_matrix(reorder_matrix(W, part), partition = part, main = "Reordered")
#'
#' @export
plot_matrix <- function(z,
                        partition = NULL,
                        zlim = NULL,
                        col = NULL,
                        max_cells = 2000L,
                        border = "grey20",
                        file = NULL,
                        format = NULL,
                        width = 2500,
                        height = 2350,
                        res = 300,
                        main = NULL,
                        xlab = "",
                        ylab = "",
                        ...) {
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z) || !is.numeric(z)) {
    abort_icons("`z` must be a numeric matrix.", class = "icons_type_error")
  }
  if (!all(is.finite(z))) {
    abort_icons("`z` must contain only finite values.", class = "icons_na_error")
  }

  sizes <- NULL
  if (!is.null(partition)) {
    if (!inherits(partition, "icons_partition")) {
      abort_icons("`partition` must be an <icons_partition>.",
                  class = "icons_type_error")
    }
    sizes <- partition$sizes
  }

  if (is.null(zlim)) {
    zlim <- range(z)
  } else if (!is.numeric(zlim) || length(zlim) != 2L) {
    abort_icons("`zlim` must be a numeric vector of length 2.",
                class = "icons_type_error")
  }
  if (is.null(col)) col <- jet_colors(256L)

  # Block-average oversized matrices so drawing stays proportional to what the
  # device can actually show.  The block factor is a whole number of cells, so
  # original index i always lands in cell ceiling(i / blk) and the axis mapping
  # below is exact.
  dim_full <- dim(z)
  blk <- 1L
  if (is.finite(max_cells) && max(dim_full) > max_cells) {
    blk <- as.integer(ceiling(max(dim_full) / max_cells))
    z <- downsample(z, blk)
  }

  if (!is.null(file)) {
    format <- if (is.null(format)) tolower(tools::file_ext(file)) else tolower(format)
    if (format == "jpg") format <- "jpeg"
    open_device(file, format, width, height, res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(op)
    graphics::layout(1)
  }, add = TRUE, after = FALSE)

  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(5.5, 1))
  graphics::par(mar = c(4.5, 4.5, if (is.null(main)) 2 else 3.5, 1))

  nr <- nrow(z)
  nc <- ncol(z)
  graphics::image(x = seq_len(nc), y = seq_len(nr),
                  z = t(z[nr:1L, , drop = FALSE]),
                  zlim = zlim, col = col, axes = FALSE, ann = FALSE, ...)

  # Ticks are chosen in the original index range and then mapped onto the drawn
  # grid.  Doing it the other way round -- pretty() on the drawn grid with the
  # labels rescaled afterwards -- puts the ticks at round *positions* but gives
  # ragged labels such as 1050, 2100, 3149.
  at_x <- index_ticks(dim_full[2L])
  graphics::axis(1, at = at_x / blk, labels = at_x)
  at_y <- index_ticks(dim_full[1L])
  graphics::axis(2, at = nr - at_y / blk + 1, labels = at_y, las = 1)
  graphics::title(main = main, xlab = xlab, ylab = ylab)

  if (!is.null(sizes) && !is.na(border) && length(sizes) > 1L) {
    b <- cumsum(sizes) / blk
    b <- b[b < min(nr, nc)]
    graphics::abline(v = b + 0.5, col = border, lwd = 0.7)
    graphics::abline(h = nr - b + 0.5, col = border, lwd = 0.7)
  }
  graphics::box()

  # Colour bar.
  graphics::par(mar = c(4.5, 0.5, if (is.null(main)) 2 else 3.5, 4))
  graphics::plot(NA, xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i",
                 ann = FALSE, axes = FALSE)
  yy <- seq(zlim[1L], zlim[2L], length.out = length(col))
  graphics::rect(0, utils::head(yy, -1L), 1, utils::tail(yy, -1L),
                 col = col[-1L], border = NA)
  graphics::axis(4, at = bar_ticks(zlim), las = 2, lwd = 0, line = -0.6)
  graphics::box()

  if (!is.null(file)) {
    message("Saved to ", normalizePath(file, mustWork = FALSE))
  }
  invisible(NULL)
}

jet_colors <- function(n) {
  grDevices::colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
      "yellow", "#FF7F00", "red", "#7F0000")
  )(n)
}

# Nicely rounded tick positions over the index range 1..n, MATLAB style: round
# numbers only, and never a tick outside the data.
index_ticks <- function(n) {
  at <- pretty(c(1, n))
  at[at >= 1 & at <= n]
}

# Colour-bar ticks.  Dropping the pretty() values that fall outside zlim can
# leave the bar sparsely labelled -- pretty(c(-0.48, 1), 5) keeps only 0, 0.5
# and 1 -- so step up the requested density until enough ticks survive.
bar_ticks <- function(zlim) {
  if (!is.finite(diff(zlim)) || diff(zlim) <= 0) return(zlim[1L])
  at <- numeric(0)
  for (k in c(5L, 6L, 8L, 10L)) {
    at <- pretty(zlim, k)
    at <- at[at >= zlim[1L] & at <= zlim[2L]]
    if (length(at) >= 4L) break
  }
  at
}

# Block-average with an integer block size, so every cell but possibly the last
# covers exactly `blk` original rows and columns.
downsample <- function(z, blk) {
  ri <- rep(seq_len(ceiling(nrow(z) / blk)), each = blk, length.out = nrow(z))
  ci <- rep(seq_len(ceiling(ncol(z) / blk)), each = blk, length.out = ncol(z))
  m <- rowsum(z, ri, reorder = TRUE)
  m <- t(rowsum(t(m), ci, reorder = TRUE))
  m / outer(tabulate(ri), tabulate(ci))
}

open_device <- function(file, format, width, height, res) {
  switch(format,
    tiff = grDevices::tiff(file, width = width, height = height, res = res,
                           compression = "lzw"),
    png  = grDevices::png(file, width = width, height = height, res = res),
    jpeg = grDevices::jpeg(file, width = width, height = height, res = res),
    svg  = grDevices::svg(file, width = width / res, height = height / res),
    eps  = grDevices::postscript(file, width = width / res, height = height / res,
                                 paper = "special", horizontal = FALSE,
                                 onefile = FALSE),
    pdf  = grDevices::pdf(file, width = width / res, height = height / res),
    abort_icons("unsupported format ", sQuote(format),
                "; use tiff, png, jpeg, pdf, eps or svg.",
                class = "icons_value_error")
  )
}

#' Plot a detected partition
#'
#' Convenience wrapper: reorders `W` by the partition and draws it.
#'
#' @param x An `"icons_partition"` object.
#' @param W The matrix the partition was detected from.
#' @param ... Passed to [plot_matrix()].
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' data(sim)
#' W <- cor(sim)
#' plot(icons_detect(W, threshold = 0.6), W)
#'
#' @export
plot.icons_partition <- function(x, W, ...) {
  if (missing(W)) {
    abort_icons("`W` is required: a partition does not store the matrix.",
                class = "icons_value_error")
  }
  plot_matrix(reorder_matrix(W, x), partition = x, ...)
}
