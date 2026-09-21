art <- c(Figure_5_plasmid_comparison = "Figure_5_plasmid_comparison.png",
  Figure_S01_plasmid_network = "Figure_S1_plasmid_network.png",
  Figure_S02_confocal_biofilms = "Figure_S2_confocal_biofilms.png",
  Figure_S03_biofilm_quantification = "Figure_S3_biofilm_quantification.png",
  Figure_S04_mouse_survival = "Figure_S4_mouse_survival.png")
for (stem in names(art)) {
  im <- magick::image_read(file.path(DATA, "artwork", art[[stem]]))
  area <- if (startsWith(stem, "Figure_5")) "main" else "supplementary"
  magick::image_write(im, file.path(OUT, area, paste0(stem, ".tiff")),
    format = "tiff", density = "600x600", compression = "lzw")
}
