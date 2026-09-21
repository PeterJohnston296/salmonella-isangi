project_root <- function() {
  p <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "salmonella-isangi.Rproj"))) return(p)
    q <- dirname(p)
    if (q == p) stop("Run from the salmonella-isangi project directory.", call. = FALSE)
    p <- q
  }
}
ROOT <- project_root()
DATA <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "results", "rebuilt", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
for (x in c("main", "supplementary", "source_data", "matrices", "audit", "previews"))
  dir.create(file.path(OUT, x), recursive = TRUE, showWarnings = FALSE)

native <- new.env(parent = globalenv())
sys.source(file.path(ROOT, "scripts", "lib", "figures.R"), envir = native)
cfg <- native$isangi_config(project_root = ROOT, input_root = DATA)
cfg$formats <- c("tiff", "png")
cfg$create_assemblies <- FALSE
cfg$save_grobs <- FALSE
cfg$matrix$transform <- "sqrt"
cfg$palette$timeline <- c("Neonatal unit sampling" = "#419D78", "River water" = "#3385B6", "Blood culture" = "#C94046", CSF = "#74528D", Other = "#9A9A9A")
d <- native$load_isangi(cfg)

save_drawing <- function(draw, stem, width, height, area = "main") {
  destination <- file.path(OUT, area, paste0(stem, ".tiff"))
  preview <- file.path(OUT, "previews", paste0(stem, ".png"))
  render <- function(path, dpi, type) {
    temporary <- paste0(path, ".partial")
    device <- if (type == "tiff") ragg::agg_tiff else ragg::agg_png
    args <- list(filename = temporary, width = width, height = height,
      units = "mm", res = dpi, background = "white")
    if (type == "tiff") args$compression <- "lzw"
    do.call(device, args)
    id <- grDevices::dev.cur()
    tryCatch({grid::grid.newpage(); draw(width, height)},
      finally = grDevices::dev.off(id))
    if (!file.exists(temporary) || file.info(temporary)$size == 0) stop("Empty figure: ", stem)
    if (file.exists(path)) unlink(path)
    if (!file.rename(temporary, path)) stop("Could not save figure: ", stem)
  }
  render(destination, 600, "tiff")
  render(preview, 160, "png")
  invisible(destination)
}
save_native <- function(f, stem, area = "supplementary") {
  save_drawing(f$draw, stem, f$width * 25.4, f$height * 25.4, area)
}
save_plot <- function(p, stem, width, height, area = "main") {
  save_drawing(function(w, h) print(p, newpage = FALSE), stem, width, height, area)
}
