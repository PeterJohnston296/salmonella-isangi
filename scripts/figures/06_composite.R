paths <- file.path(OUT, "main", c("Figure_3A_CC25.tiff", "Figure_3B_AMR_gene_burden_600dpi.tiff"))
images <- lapply(paths, function(p) magick::image_trim(magick::image_background(magick::image_read(p), "white", flatten = TRUE), fuzz = 0))
widths <- c(4000, 3314)
canvas_width <- 4240L
margin <- 70L
panels <- lapply(seq_along(images), function(i) {
  im <- magick::image_resize(images[[i]], paste0(widths[i], "x"))
  im <- magick::image_border(im, "white", paste0(margin, "x", margin))
  info <- magick::image_info(im)
  stopifnot(info$width <= canvas_width)
  im <- magick::image_extent(im, paste0(canvas_width, "x", info$height + 100), gravity = "south", color = "white")
  magick::image_annotate(im, c("A", "B")[i], gravity = "northwest", location = "+35+20", size = 72, weight = 700, font = if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "DejaVu-Sans", color = "black")
})
figure <- magick::image_append(c(panels[[1]], magick::image_blank(canvas_width, 60, "white"), panels[[2]]), stack = TRUE)
magick::image_write(figure, file.path(OUT, "main", "Figure_3A_B_Genome_Medicine_600dpi.tiff"),
  format = "tiff", density = "600x600", compression = "lzw")
magick::image_write(magick::image_resize(figure, "1600x"), file.path(OUT, "previews", "Figure_3A_B_preview.png"))
