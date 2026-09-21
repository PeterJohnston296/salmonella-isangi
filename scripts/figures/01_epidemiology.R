sys.source(file.path(ROOT, "scripts", "lib", "figure1.R"), envir = environment())
r <- build_figure1(project_root = ROOT, output_dir = file.path(OUT, "figure1"), map_mode = "artwork", open_result = FALSE)
file.copy(r$tiff, file.path(OUT, "main", "Figure_1A_D_Genome_Medicine_600dpi.tiff"), overwrite = TRUE)
file.copy(r$preview, file.path(OUT, "previews", "Figure_1A_D_preview.png"), overwrite = TRUE)
