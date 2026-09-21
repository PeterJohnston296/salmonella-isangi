build_figure1 <- function(
  project_root = getwd(),
  output_dir = file.path(project_root, "results", "rebuilt", "figure1"),
  map_mode = "artwork",
  boundary_source_crs = NULL,
  river_source_crs = NULL,
  missing_geographic_crs = 4326,
  missing_projected_boundary_crs = 20936,
  missing_projected_river_crs = 21036,
  open_result = FALSE
) {

  pkgs <- c("ggplot2", "cowplot", "igraph", "ragg", "magick")
  if (map_mode == "geography") pkgs <- c(pkgs, "sf")
  absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Packages unavailable: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages(library(ggplot2))

  if (!dir.exists(project_root)) stop("Project folder not found: ", project_root)
  project_root <- normalizePath(project_root, mustWork = TRUE)
  stamp <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  out <- output_dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out)) stop("Cannot create output directory: ", out)
  audit <- character()
  record <- function(...) {
    s <- paste0(...)
    message(s)
    audit <<- c(audit, s)
  }
  on.exit(writeLines(audit, file.path(out, "build_log.txt")), add = TRUE)
  if (map_mode == "geography") {
    old_axis <- sf::st_axis_order()
    sf::st_axis_order(FALSE)
    on.exit(sf::st_axis_order(old_axis), add = TRUE)
  }

  inventory <- list.files(file.path(project_root, "data"), recursive = TRUE, full.names = TRUE)
  inventory <- inventory[!grepl("__MACOSX|[.]Rproj[.]user|figure1_submission_", inventory)]
  bare_name <- function(x) {
    ext <- tolower(tools::file_ext(x))
    stem <- tools::file_path_sans_ext(basename(x))
    stem <- sub("( *\\([0-9]+\\))+$", "", stem, perl = TRUE)
    paste0(tolower(stem), ".", ext)
  }
  find_file <- function(name, required = TRUE) {
    hits <- inventory[bare_name(inventory) == bare_name(name)]
    if (!length(hits)) {
      if (required) stop("Missing input: ", name, "\nSearch root: ", project_root)
      return(NA_character_)
    }

    rank <- ifelse(grepl("/inputs/|/figure_inputs/", hits), 0L,
                   ifelse(grepl("/outputs/", hits), 2L, 1L))
    hits <- hits[order(rank, -as.numeric(file.info(hits)$mtime), hits)]
    record(name, " -> ", hits[1])
    hits[1]
  }
  read_table <- function(path) {
    as.data.frame(utils::read.csv(path, check.names = FALSE,
      stringsAsFactors = FALSE, na.strings = c("", "NA")))
  }
  clean_name <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
  column <- function(df, alternatives, required = TRUE) {
    for (a in alternatives) {
      j <- which(clean_name(names(df)) == clean_name(a))
      if (length(j) == 1L) return(df[, j, drop = TRUE])
      if (length(j) > 1L) stop("Ambiguous column: ", a)
    }
    if (required) stop("Missing column: ", paste(alternatives, collapse = " / "),
      "\nAvailable: ", paste(names(df), collapse = ", "))
    NULL
  }
  num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
  month_date <- function(x) {
    s <- trimws(as.character(x))
    s[grepl("^[0-9]{4}-[0-9]{2}$", s)] <-
      paste0(s[grepl("^[0-9]{4}-[0-9]{2}$", s)], "-01")
    d <- as.Date(s, format = "%Y-%m-%d")
    k <- is.na(d) & !is.na(s)
    d[k] <- as.Date(s[k], format = "%d/%m/%Y")
    if (anyNA(d)) stop("Unparseable timeline dates: ", paste(unique(s[is.na(d)]), collapse = ", "))
    as.Date(format(d, "%Y-%m-01"), format = "%Y-%m-%d")
  }
  source_name <- function(x) {
    s <- tolower(trimws(as.character(x))); s[is.na(s)] <- ""
    z <- rep("Other", length(s))
    z[grepl("neonatal|chatinkha|nursery", s)] <- "Neonatal unit sampling"
    z[grepl("river", s)] <- "River water"
    z[grepl("blood", s)] <- "Blood culture"
    z[grepl("csf|cerebrospinal", s)] <- "CSF"
    z
  }
  source_cols <- c("Neonatal unit sampling" = "#4DAA7C", "River water" = "#4E93C6",
                   "Blood culture" = "#C84D5A", "CSF" = "#775098", "Other" = "#858B92")
  font <- "sans"
  pub_theme <- function() {
    theme_classic(base_size = 7.5, base_family = font) + theme(
      axis.text = element_text(colour = "#30353B", size = 6.8),
      axis.title = element_text(size = 7.5),
      axis.line = element_line(linewidth = 0.28, colour = "#30353B"),
      axis.ticks = element_line(linewidth = 0.25),
      axis.ticks.length = grid::unit(1, "mm"),
      legend.position = "none", plot.margin = margin(2, 2, 2, 2),
      plot.background = element_rect(fill = "white", colour = NA))
  }

  annual_path <- find_file("BSI_annual_resistance_summary_2015_2023.csv")
  time_path <- find_file("Malawi_timeline_74_month_resolution.csv")
  time_raw <- read_table(time_path)
  timeline <- data.frame(
    isolate_id = trimws(as.character(column(time_raw, "isolate_id"))),
    source = source_name(column(time_raw, "source")),
    collection_month = month_date(column(time_raw, "collection_month")),
    stringsAsFactors = FALSE)
  if (nrow(timeline) != 74L || anyDuplicated(timeline$isolate_id))
    stop("Expected 74 unique Malawi timeline isolates.")

  matrix_names <- c("CC25_complete_331.csv", "CC25_tree_aligned_327.csv",
                    "ST335_tree_aligned_224.csv", "Malawi_74_source_matrix.csv")
  D <- NULL; matrix_path <- NULL
  for (nm in matrix_names) {
    p <- find_file(nm, required = FALSE)
    if (is.na(p)) next
    candidate <- tryCatch({
      tab <- read_table(p)
      ids <- trimws(as.character(tab[, 1, drop = TRUE]))
      m <- as.matrix(tab[, -1, drop = FALSE])
      storage.mode(m) <- "numeric"
      rownames(m) <- ids; colnames(m) <- trimws(colnames(m))
      wanted <- timeline$isolate_id
      if (anyDuplicated(ids) || anyDuplicated(colnames(m)) ||
          !all(wanted %in% ids) || !all(wanted %in% colnames(m))) stop("ID mismatch")
      m <- m[wanted, wanted, drop = FALSE]
      if (any(!is.finite(m)) || any(m < 0) || any(m != round(m)) ||
          any(diag(m) != 0) || any(m != t(m))) stop("Invalid SNP matrix")
      m
    }, error = function(e) { record("Rejected matrix ", p, ": ", conditionMessage(e)); NULL })
    if (!is.null(candidate)) { D <- candidate; matrix_path <- p; break }
  }
  if (is.null(D)) stop("No valid, symmetric SNP matrix covering all 74 isolates was found.")
  reference <- "CAAP2A"
  if (!reference %in% colnames(D)) stop("Reference CAAP2A is absent from the matrix.")
  timeline$SNP <- as.numeric(D[timeline$isolate_id, reference])
  record("Validated Malawi matrix: ", nrow(D), " isolates; range 0-", max(D), " SNPs.")
  old_dist <- column(time_raw, "SNP_distance_to_reference", required = FALSE)
  if (!is.null(old_dist) && any(num(old_dist) != timeline$SNP, na.rm = TRUE))
    record("Timeline distances differ from the matrix; the validated matrix supplies plotted distances.")

  ar <- read_table(annual_path)
  a <- data.frame(year = num(column(ar, "year")), group = toupper(trimws(column(ar, "resistance_group"))),
                  n = num(column(ar, "n")), total = num(column(ar, "annual_total")))
  if (any(!is.finite(a$year)) || any(!is.finite(a$n)) || any(a$n < 0) ||
      any(!is.finite(a$total)) || any(a$total <= 0)) stop("Invalid annual resistance table.")
  a_total <- unique(a[, c("year", "total")])
  if (anyDuplicated(a_total$year)) stop("Conflicting annual totals.")
  a_total <- a_total[order(a_total$year), ]
  a_total$MDR <- vapply(a_total$year, function(y) sum(a$n[a$year == y & a$group == "MDR"]), numeric(1))
  a_total$XDR <- vapply(a_total$year, function(y) sum(a$n[a$year == y & a$group == "XDR"]), numeric(1))
  totals_check <- vapply(a_total$year, function(y) sum(a$n[a$year == y]), numeric(1))
  if (any(totals_check != a_total$total)) stop("Annual group counts do not sum to annual totals.")
  a_total$p_mdr <- a_total$MDR / a_total$total
  a_total$p_xdr <- a_total$XDR / a_total$total
  count_top <- ceiling(max(a_total$total) / 20) * 20
  prop_top <- max(0.6, ceiling(max(a_total$p_mdr, a_total$p_xdr) * 10) / 10)
  factor_a <- count_top / prop_top
  pA <- ggplot(a_total, aes(x = year)) +
    geom_col(aes(y = total), width = 0.68, fill = "#E0E5E9", colour = "#98A3AE", linewidth = 0.25) +
    geom_line(aes(y = p_mdr * factor_a), colour = "#DB9345", linewidth = 0.55) +
    geom_point(aes(y = p_mdr * factor_a), colour = "#DB9345", size = 1.55) +
    geom_line(aes(y = p_xdr * factor_a), colour = "#BD4B68", linewidth = 0.55) +
    geom_point(aes(y = p_xdr * factor_a), colour = "#BD4B68", size = 1.55) +
    scale_x_continuous(breaks = a_total$year, expand = expansion(add = 0.55)) +
    scale_y_continuous(limits = c(0, count_top), breaks = seq(0, count_top, 40),
      expand = expansion(mult = c(0, 0.02)),
      sec.axis = sec_axis(~ . / factor_a, name = "Isolates (%)",
        breaks = seq(0, prop_top, 0.2), labels = function(x) paste0(round(100 * x), "%"))) +
    labs(x = "Collection year", y = "NTS blood-culture isolates") + pub_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y.right = element_text(margin = margin(l = 3)))

  bubbles <- aggregate(list(n = rep(1L, nrow(timeline))),
    timeline[, c("collection_month", "source", "SNP")], sum)

  bubbles$plot_date <- bubbles$collection_month + 14
  max_bubble <- max(bubbles$n)
  timeline_plot <- function(z, detailed = FALSE) {
    breaks <- if (detailed) as.Date(paste0("2020-", c("01", "03", "05", "07", "09", "11"), "-15")) else
      as.Date(paste0(2018:2023, "-07-01"))
    xl <- if (detailed) as.Date(c("2020-01-01", "2020-12-31")) else
      as.Date(c("2018-01-01", "2023-08-31"))
    ymax <- if (detailed) max(9, max(z$SNP) + 1) else ceiling(max(z$SNP) / 5) * 5
    ggplot(z, aes(plot_date, SNP)) +
      geom_point(aes(fill = source, size = n), shape = 21, colour = "white", stroke = 0.18, alpha = 0.94) +
      scale_fill_manual(values = source_cols) +
      scale_size_area(max_size = 5.2, limits = c(0, max_bubble)) +
      scale_x_date(breaks = breaks, date_labels = if (detailed) "%b" else "%Y",
                   limits = xl, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(limits = c(0, ymax), breaks = if (detailed) seq(0, ymax, 2) else seq(0, ymax, 5),
                         expand = expansion(mult = c(0.03, 0.06))) +
      labs(x = if (detailed) "Collection month (2020)" else NULL, y = NULL) + pub_theme()
  }
  pBi <- timeline_plot(bubbles)
  b2020 <- bubbles[format(bubbles$collection_month, "%Y") == "2020", ]
  if (!nrow(b2020)) stop("No 2020 observations.")
  pBii <- timeline_plot(b2020, TRUE)

  group <- rep(NA_integer_, nrow(D)); reps <- integer(); ng <- 0L
  for (i in seq_len(nrow(D))) if (is.na(group[i])) {
    same <- which(vapply(seq_len(nrow(D)), function(j) all(D[i, ] == D[j, ]), logical(1)))
    if (any(D[same, same, drop = FALSE] != 0)) stop("Invalid zero-distance grouping.")
    ng <- ng + 1L; group[same] <- ng; reps <- c(reps, i)
  }
  if (ng < 2) stop("Network needs at least two distinct profiles.")
  Dr <- D[reps, reps, drop = FALSE]
  ix <- which(upper.tri(Dr), arr.ind = TRUE)
  ord <- order(Dr[ix], ix[, 1], ix[, 2])
  parent <- seq_len(ng); ed <- matrix(NA_real_, ng - 1L, 3); k <- 0L
  root_of <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (z in ord) {
    u <- ix[z, 1]; v <- ix[z, 2]; ru <- root_of(u); rv <- root_of(v)
    if (ru != rv) {
      parent[rv] <- ru; k <- k + 1L; ed[k, ] <- c(u, v, Dr[u, v])
      if (k == ng - 1L) break
    }
  }
  if (k != ng - 1L) stop("MST is not connected.")
  edges <- as.data.frame(ed); names(edges) <- c("from", "to", "SNP")
  graph <- igraph::graph_from_data_frame(
    data.frame(from = as.character(edges$from), to = as.character(edges$to)),
    directed = FALSE, vertices = data.frame(name = as.character(seq_len(ng))))
  xy <- igraph::layout_with_kk(graph, coords = igraph::layout_in_circle(graph),
                              weights = rep(1, nrow(edges)), maxiter = 5000)

  xy <- stats::prcomp(xy, center = TRUE, scale. = FALSE)$x
  Wc <- 80; Hc <- 73
  for (j in 1:2) xy[, j] <- (xy[, j] - min(xy[, j])) / max(1e-8, diff(range(xy[, j])))
  xy[, 1] <- 7 + xy[, 1] * (Wc - 14)
  xy[, 2] <- 7 + xy[, 2] * (Hc - 14)
  node_n <- tabulate(group, nbins = ng)
  radius <- 0.76 * sqrt(node_n)

  anchor <- xy
  for (iteration in 1:200) {
    shift <- matrix(0, ng, 2)
    for (i in 1:(ng - 1L)) for (j in (i + 1L):ng) {
      delta <- xy[j, ] - xy[i, ]; len <- sqrt(sum(delta^2))
      minimum <- radius[i] + radius[j] + 5.5
      if (len < minimum) {
        direction <- if (len > 1e-8) delta / len else c(cos(i + j), sin(i + j))
        push <- direction * (minimum - len) * 0.24
        shift[i, ] <- shift[i, ] - push; shift[j, ] <- shift[j, ] + push
      }
    }
    xy <- xy + shift + 0.004 * (anchor - xy)
    xy[, 1] <- pmax(5, pmin(Wc - 5, xy[, 1]))
    xy[, 2] <- pmax(5, pmin(Hc - 5, xy[, 2]))
  }
  member <- data.frame(group_id = sprintf("G%02d", group), isolate_id = timeline$isolate_id,
                        source = timeline$source, stringsAsFactors = FALSE)
  nodes <- data.frame(group_id = sprintf("G%02d", seq_len(ng)), n = node_n,
                       x_mm = xy[, 1], y_mm = xy[, 2])
  record("MST: ", ng, " display groups, ", nrow(edges), " edges; zero-SNP edges retained.")

  grobs <- list(); add <- function(g) { grobs[length(grobs) + 1L] <<- list(g) }
  ux <- function(x) grid::unit(x / Wc, "npc")
  uy <- function(y) grid::unit(y / Hc, "npc")
  for (i in seq_len(nrow(edges))) {
    u <- edges$from[i]; v <- edges$to[i]
    add(grid::segmentsGrob(ux(xy[u, 1]), uy(xy[u, 2]), ux(xy[v, 1]), uy(xy[v, 2]),
                          gp = grid::gpar(col = "#687680", lwd = 0.8)))
  }
  for (i in seq_len(ng)) {
    counts <- table(factor(timeline$source[group == i], levels = names(source_cols)))
    start <- 0
    for (j in seq_along(counts)) if (counts[j] > 0) {
      end <- start + 2 * pi * as.numeric(counts[j]) / sum(counts)
      angle <- seq(start, end, length.out = max(12, ceiling(80 * (end - start) / (2 * pi))))
      add(grid::polygonGrob(ux(c(xy[i, 1], xy[i, 1] + radius[i] * cos(angle))),
        uy(c(xy[i, 2], xy[i, 2] + radius[i] * sin(angle))),
        gp = grid::gpar(fill = unname(source_cols[j]), col = NA)))
      start <- end
    }
    add(grid::circleGrob(ux(xy[i, 1]), uy(xy[i, 2]), r = grid::unit(radius[i], "mm"),
                        gp = grid::gpar(fill = NA, col = "#45515D", lwd = 0.55)))
  }

  boxes <- data.frame(x = numeric(), y = numeric(), w = numeric(), h = numeric())
  place_label <- function(candidates, label, fontsize, colour, bold = FALSE) {
    bw <- nchar(label) * fontsize * 0.19 + 0.8; bh <- fontsize * 0.353 + 0.65
    score <- apply(candidates, 1, function(p) {
      dx <- pmax(abs(xy[, 1] - p[1]) - bw / 2, 0)
      dy <- pmax(abs(xy[, 2] - p[2]) - bh / 2, 0)
      hit_nodes <- sum(dx^2 + dy^2 < (radius + 0.35)^2)
      hit_labels <- if (nrow(boxes)) sum(abs(boxes$x - p[1]) < (boxes$w + bw) / 2 &
                                          abs(boxes$y - p[2]) < (boxes$h + bh) / 2) else 0
      outside <- (p[1] - bw / 2 < 0 || p[1] + bw / 2 > Wc ||
                  p[2] - bh / 2 < 0 || p[2] + bh / 2 > Hc)
      1000 * outside + 100 * hit_nodes + 100 * hit_labels
    }) + seq_len(nrow(candidates)) * 0.001
    p <- candidates[which.min(score), ]
    boxes <<- rbind(boxes, data.frame(x = p[1], y = p[2], w = bw, h = bh))
    add(grid::roundrectGrob(ux(p[1]), uy(p[2]), width = grid::unit(bw, "mm"),
      height = grid::unit(bh, "mm"), r = grid::unit(0.25, "mm"),
      gp = grid::gpar(fill = "#FFFFFFED", col = NA)))
    add(grid::textGrob(label, ux(p[1]), uy(p[2]), gp = grid::gpar(
      fontsize = fontsize, fontfamily = font, fontface = if (bold) "bold" else "plain", col = colour)))
  }
  for (i in seq_len(nrow(edges))) {
    u <- edges$from[i]; v <- edges$to[i]; dv <- xy[v, ] - xy[u, ]
    len <- max(1e-8, sqrt(sum(dv^2))); normal <- c(-dv[2], dv[1]) / len
    start <- xy[u, ] + dv / len * (radius[u] + 0.3)
    end <- xy[v, ] - dv / len * (radius[v] + 0.3)
    candidates <- do.call(rbind, lapply(c(0, 1.2, -1.2, 2.4, -2.4, 3.6, -3.6), function(off)
      t(vapply(c(0.5, 0.3, 0.7, 0.15, 0.85), function(f) start + f * (end - start) + off * normal, numeric(2)))))
    place_label(candidates, as.character(edges$SNP[i]), 5.8, "#25313C", TRUE)
  }
  network_grob <- do.call(grid::grobTree, grobs)

  bp <- rp <- ep <- ip <- NA_character_
  other_points <- isangi_points <- data.frame()
  if (map_mode == "geography") {

  map_dir <- file.path(project_root, "map_data")
  bp <- file.path(map_dir, "dist_bnd.shp"); rp <- file.path(map_dir, "rivers2.shp")
  ep <- file.path(map_dir, "ERST_summary.csv"); ip <- file.path(map_dir, "isangi_metadata.csv")
  must_exist <- c(bp, rp, ep, ip)
  if (any(!file.exists(must_exist))) stop("Missing map files:\n", paste(must_exist[!file.exists(must_exist)], collapse = "\n"))
  map_crs <- 32736
  read_spatial <- function(path, override, projected_default) {
    x <- sf::st_read(path, quiet = TRUE)
    if (!nrow(x)) stop("Empty spatial file: ", path)
    b <- sf::st_bbox(x)
    if (any(!is.finite(b))) stop("Invalid bounding box: ", path)
    geographic <- b["xmin"] >= 30 && b["xmax"] <= 38 && b["ymin"] >= -20 && b["ymax"] <= -7
    record(basename(path), " raw bounds: ", paste(round(b, 5), collapse = ", "))
    if (!is.null(override)) {
      sf::st_crs(x) <- override
      record("User-specified source CRS: ", override)
    } else if (geographic && (is.na(sf::st_crs(x)) || !isTRUE(sf::st_is_longlat(x)))) {
      sf::st_crs(x) <- missing_geographic_crs
      record("ASSUMPTION: ", basename(path), " contains longitude/latitude; source CRS set to ",
             missing_geographic_crs, ". Its datum is not proved by the coordinate range.")
    } else if (is.na(sf::st_crs(x))) {
      if (b["xmin"] < 100000 || b["xmax"] > 1100000 || b["ymin"] < 7000000 || b["ymax"] > 9300000)
        stop("Unrecognised coordinate units in ", basename(path), ". Supply its true source CRS.")
      sf::st_crs(x) <- projected_default
      record("ASSUMPTION: projected source CRS for ", basename(path), " = ", projected_default)
    }
    x <- sf::st_transform(x, map_crs)
    x <- sf::st_make_valid(x)
    x <- x[!sf::st_is_empty(x), ]
    if (!nrow(x)) stop("No geometries survived transformation: ", path)
    x
  }
  point_table <- function(path, is_isangi = FALSE) {
    raw <- read_table(path)
    z <- data.frame(longitude = num(column(raw, c("Longitude", "Lon", "Long"))),
                    latitude = num(column(raw, c("Latitude", "Lat"))))
    valid <- is.finite(z$longitude) & is.finite(z$latitude)
    if (is_isangi && any(!valid)) stop("Missing Isangi map coordinates: resolve before submission.")
    if (is_isangi) keep <- valid else {
      w <- toupper(trimws(as.character(column(raw, c("WGS Result", "Serovar")))))
      keep <- valid & !is.na(w) & !w %in% c("", "NT", "NA", "N/A", "NOT TESTED", "UNKNOWN") &
        !grepl("ISANGI", w) & !grepl("^(SALMONELLA[ _]+)?TYPHI$", w)
    }
    z <- z[keep, , drop = FALSE]
    if (!nrow(z)) stop("No usable sampling coordinates in ", basename(path))
    if (any(z$longitude < 34 | z$longitude > 36 | z$latitude < -17 | z$latitude > -15))
      stop("Coordinates outside the expected Blantyre region in ", basename(path), ". No points were silently discarded.")
    z$source_row <- which(keep)
    record(basename(path), ": ", nrow(z), " plotted records; ", nrow(unique(z[, 1:2])), " distinct coordinate pairs.")
    z
  }
  other_points <- point_table(ep)
  isangi_points <- point_table(ip, TRUE)
  other_sf <- sf::st_transform(sf::st_as_sf(other_points, coords = c("longitude", "latitude"), crs = 4326), map_crs)
  isangi_sf <- sf::st_transform(sf::st_as_sf(isangi_points, coords = c("longitude", "latitude"), crs = 4326), map_crs)
  qech <- sf::st_transform(sf::st_as_sf(data.frame(longitude = 35.021550, latitude = -15.803036),
                             coords = c("longitude", "latitude"), crs = 4326), map_crs)
  boundary <- read_spatial(bp, boundary_source_crs, missing_projected_boundary_crs)
  rivers <- read_spatial(rp, river_source_crs, missing_projected_river_crs)

  all_xy <- rbind(sf::st_coordinates(other_sf), sf::st_coordinates(isangi_sf), sf::st_coordinates(qech))
  xr <- range(all_xy[, 1]); yr <- range(all_xy[, 2])
  ww <- max(diff(xr) * 1.18, 10000); hh <- max(diff(yr) * 1.18, 11000)
  frame_ratio <- 80 / 73
  if (ww / hh < frame_ratio) ww <- hh * frame_ratio else hh <- ww / frame_ratio
  if (ww > 70000 || hh > 70000 || ww < 1000 || hh < 1000) stop("Implausible map dimensions.")
  xc <- mean(xr); yc <- mean(yr)
  bb <- sf::st_bbox(c(xmin = xc - ww/2, ymin = yc - hh/2, xmax = xc + ww/2, ymax = yc + hh/2), crs = sf::st_crs(map_crs))
  boundary <- suppressWarnings(sf::st_crop(boundary, bb))
  rivers <- suppressWarnings(sf::st_crop(rivers, bb))
  boundary <- boundary[!sf::st_is_empty(boundary), ]
  rivers <- rivers[!sf::st_is_empty(rivers), ]
  if (!nrow(boundary) || !nrow(rivers))
    stop("Map layers do not overlap the sampling region. Check source CRS; an empty map will not be exported.")
  land <- sf::st_union(sf::st_geometry(boundary))
  distance_to_land <- min(as.numeric(sf::st_distance(qech, land)))
  if (!is.finite(distance_to_land) || distance_to_land > 2000)
    stop("QECH is not within/near the boundary: source CRS or boundary selection is wrong.")
  river_length <- sum(as.numeric(sf::st_length(rivers)), na.rm = TRUE)
  if (river_length < 500) stop("Less than 500 m of waterways in the map. Check source CRS.")
  record("Map frame: ", round(ww/1000, 1), " x ", round(hh/1000, 1), " km; ", nrow(rivers), " waterway features.")
  flow <- column(sf::st_drop_geometry(rivers), c("FLOW_M3", "FLOWM3"), required = FALSE)
  rivers$width_mm <- if (is.null(flow)) rep(0.19, nrow(rivers)) else {
    f <- num(flow); ifelse(!is.finite(f), 0.14, ifelse(f > 2, 0.40, ifelse(f > 0.5, 0.28, 0.16)))
  }

  names_r <- column(sf::st_drop_geometry(rivers), c("RIVER_NAME", "RIVER_NAM", "RIV_NAME", "NAME", "RIVER"), required = FALSE)
  river_labels <- NULL
  if (!is.null(names_r)) {
    names_r <- trimws(as.character(names_r))
    good <- !is.na(names_r) & grepl("[A-Za-z]", names_r) & !tolower(names_r) %in% c("unknown", "unnamed", "river", "other")
    rr <- rivers[good, ]; rr$waterway_name <- names_r[good]
    if (nrow(rr)) {
      rr <- rr[order(as.numeric(sf::st_length(rr)), decreasing = TRUE), ]
      rr <- rr[!duplicated(tolower(rr$waterway_name)), ]
      rr <- head(rr, 3)
      river_labels <- suppressWarnings(sf::st_point_on_surface(rr))
      if (!is.null(river_labels) && nrow(river_labels) > 0 && nrow(isangi_sf) > 0) {
        dmat <- sf::st_distance(river_labels, isangi_sf)
        min_dist <- apply(dmat, 1, min, na.rm = TRUE)
        river_labels <- river_labels[as.numeric(min_dist) > 450, ]
      }
      if (!is.null(river_labels) && nrow(river_labels) > 0) {
        rxy <- sf::st_coordinates(river_labels)
        river_labels$x <- rxy[, 1]
        river_labels$y <- rxy[, 2]
        river_labels$nudge_x <- rep(0, nrow(river_labels))
        river_labels$nudge_y <- c(220, -220, 180)[seq_len(nrow(river_labels))]
      }
    }
  }
  if (is.null(river_labels) || nrow(river_labels) == 0) record("No verified river-name attributes safely placeable away from Isangi points: no river names drawn.")
  qxy <- sf::st_coordinates(qech)[1, ]
  qlab <- data.frame(x = qxy[1] + ww * 0.10, y = qxy[2] + hh * 0.075, label = "QECH")
  land_fill <- "#F0F4ED"
  pD <- ggplot() +
    annotate("rect", xmin = bb["xmin"], xmax = bb["xmax"], ymin = bb["ymin"], ymax = bb["ymax"],
             fill = land_fill, colour = NA) +
    geom_sf(data = boundary, fill = land_fill, colour = "#BAC5B3", linewidth = 0.22) +
    geom_sf(data = rivers, aes(linewidth = width_mm), colour = "#75AFC5") + scale_linewidth_identity() +
    geom_sf(data = other_sf, shape = 16, size = 1.05, colour = "#8A949B", alpha = 0.55) +
    geom_sf(data = isangi_sf, shape = 23, size = 1.9, fill = "#B63978", colour = "white", stroke = 0.18) +
    geom_sf(data = qech, shape = 23, size = 2.6, fill = "#E4A340", colour = "#2D3339", stroke = 0.35) +
    annotate("segment", x = qxy[1], y = qxy[2], xend = qlab$x, yend = qlab$y,
             linewidth = 0.25, colour = "#505C65") +
    geom_label(data = qlab, aes(x, y, label = label), size = 2.5, fontface = "bold",
               fill = "white", colour = "#25313B", label.padding = grid::unit(0.5, "mm"), linewidth = 0) +
    coord_sf(crs = sf::st_crs(map_crs), default_crs = sf::st_crs(map_crs),
      xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), datum = NA, expand = FALSE) +
    theme_void(base_family = font) + theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = land_fill, colour = NA),
      plot.background = element_rect(fill = land_fill, colour = NA)
    )
  if (!is.null(river_labels) && nrow(river_labels) > 0) {
    river_df <- as.data.frame(river_labels)
    pD <- pD + geom_label(
      data = river_df,
      aes(x = x + nudge_x, y = y + nudge_y, label = waterway_name),
      size = 2.15, fontface = "italic", family = font,
      colour = "#376C82", fill = scales::alpha(land_fill, 0.92),
      label.size = 0, label.padding = grid::unit(0.10, "lines")
    )
  }
  scale_m <- max(c(500, 1000, 2000, 5000, 10000)[c(500, 1000, 2000, 5000, 10000) <= ww * 0.24])
  sx <- unname(bb["xmin"]) + ww * 0.055; sy <- unname(bb["ymin"]) + hh * 0.05
  pD <- pD +
    annotate("rect", xmin = sx - ww*0.012, xmax = sx + scale_m + ww*0.035,
             ymin = sy - hh*0.013, ymax = sy + hh*0.065, fill = "white", colour = NA, alpha = 0.92) +
    annotate("segment", x = sx, xend = sx + scale_m, y = sy, yend = sy, linewidth = 0.65, colour = "#26313A") +
    annotate("segment", x = c(sx, sx+scale_m), xend = c(sx, sx+scale_m),
             y = sy-hh*0.006, yend = sy+hh*0.006, linewidth = 0.5) +
    annotate("text", x = sx + scale_m/2, y = sy + hh*0.035, label = paste0(scale_m/1000, " km"),
             size = 2.3, colour = "#26313A")

  north <- sf::st_transform(sf::st_as_sf(data.frame(lon=35.021550, lat=-15.802036),
                                        coords=c("lon", "lat"), crs=4326), map_crs)
  nv <- sf::st_coordinates(north)[1, ] - qxy; nv <- nv / sqrt(sum(nv^2))
  nx <- unname(bb["xmax"]) - ww*0.07; ny <- unname(bb["ymax"]) - hh*0.13
  pD <- pD + annotate("segment", x=nx, y=ny, xend=nx+nv[1]*hh*0.055,
       yend=ny+nv[2]*hh*0.055, arrow=grid::arrow(length=grid::unit(1.6,"mm")), linewidth=0.35) +
    annotate("text", x=nx, y=ny+hh*0.078, label="N", size=2.4, colour="#26313A")

  } else {
    image <- magick::image_read(file.path(project_root, "data", "artwork", "Figure_1D_Blantyre_map.png"))
    pD <- cowplot::ggdraw() + cowplot::draw_image(image, x = 0, y = 0, width = 1, height = 1)
    record("Panel D uses the archived map artwork.")
  }

  key_grob <- function(labels, colours, shapes, xx, yy, width, height, size=7) {
    entries <- list()
    for (i in seq_along(labels)) {
      if (shapes[i] == -1) {
        symbol <- grid::segmentsGrob(grid::unit(xx[i]/width,"npc"), grid::unit(yy[i]/height,"npc"),
                    grid::unit((xx[i]+3)/width,"npc"), grid::unit(yy[i]/height,"npc"),
                    gp=grid::gpar(col=colours[i], lwd=1.6))
      } else {
        symbol <- grid::pointsGrob(grid::unit((xx[i]+1.5)/width,"npc"), grid::unit(yy[i]/height,"npc"),
                    pch=shapes[i], size=grid::unit(2.2,"mm"),
                    gp=grid::gpar(col=colours[i], fill=colours[i]))
      }
      entries[length(entries)+1L] <- list(symbol)
      entries[length(entries)+1L] <- list(grid::textGrob(labels[i],
         x=grid::unit((xx[i]+5)/width,"npc"), y=grid::unit(yy[i]/height,"npc"), just="left",
         gp=grid::gpar(fontsize=size, fontfamily=font, col="#30353B")))
    }
    do.call(grid::grobTree, entries)
  }
  akey <- key_grob(c("MDR, not XDR", "XDR"), c("#DB9345", "#BD4B68"), c(-1,-1),
                   c(0,37), c(3,3), 64,6,7.0)
  present <- names(source_cols)[names(source_cols) %in% timeline$source]
  sx_key <- if(length(present)<=4) c(0,59,96,137)[seq_along(present)] else c(0,52,84,119,141)
  skey <- key_grob(present, unname(source_cols[present]), rep(16,length(present)),
                   sx_key, rep(4,length(present)), 164,8,7.0)
  mkey <- key_grob(c("Other NTS", "S. Isangi", "QECH", "Waterways"),
          c("#8A949B", "#B63978", "#E4A340", "#75AFC5"), c(16,23,23,-1),
          c(0,29,56,0), c(8,8,8,2),80,11,6.8)

  FW <- 170; FH <- 190
  canvas <- cowplot::ggdraw()
  put_plot <- function(p,x,y,w,h) cowplot::draw_plot(p,x/FW,y/FH,w/FW,h/FH)
  put_grob <- function(p,x,y,w,h) cowplot::draw_grob(p,x/FW,y/FH,w/FW,h/FH)
  put_text <- function(s,x,y,size=7.5,face="plain",angle=0,hjust=0,vjust=0.5)
    cowplot::draw_label(s,x/FW,y/FH,size=size,fontfamily=font,fontface=face,
                       angle=angle,hjust=hjust,vjust=vjust,colour="#26313A")
  figure <- canvas +
    put_plot(pA,1,105,78,71) + put_grob(akey,14,178,64,6) +
    put_plot(pBi,94,148,75,33) + put_plot(pBii,94,105,75,34) +
    put_text("SNPs from CAAP2A",87,145,7.5,angle=90,hjust=0.5) +
    put_text("i",91,178,7.5,"bold") + put_text("ii",91,136,7.5,"bold") +
    put_grob(skey,3,94,164,8) +
    put_grob(network_grob,3,16,Wc,Hc) + put_plot(pD,88,16,80,73) +
    {if (map_mode == "geography") put_grob(mkey,88,2,80,11) else cowplot::draw_label("")} +
    put_text("Node area: isolate count; bold edge labels: pairwise SNP distance.",3,8.0,6.3) +
    put_text("A",2,188,11,"bold") + put_text("B",88,188,11,"bold") +
    put_text("C",2,92,11,"bold") + put_text("D",88,92,11,"bold")

  utils::write.csv(a_total,file.path(out,"Figure_1A_source.csv"),row.names=FALSE)
  utils::write.csv(timeline,file.path(out,"Figure_1B_isolate_data.csv"),row.names=FALSE)
  utils::write.csv(bubbles,file.path(out,"Figure_1B_bubble_data.csv"),row.names=FALSE)
  utils::write.csv(member,file.path(out,"Figure_1C_node_membership.csv"),row.names=FALSE)
  utils::write.csv(nodes,file.path(out,"Figure_1C_node_layout.csv"),row.names=FALSE)
  edge_export <- data.frame(from=sprintf("G%02d",edges$from),to=sprintf("G%02d",edges$to),SNP=edges$SNP)
  utils::write.csv(edge_export,file.path(out,"Figure_1C_edges.csv"),row.names=FALSE)
  utils::write.csv(D,file.path(out,"Malawi_74_validated_SNP_matrix.csv"))
  utils::write.csv(other_points,file.path(out,"Figure_1D_other_NTS_positions.csv"),row.names=FALSE)
  utils::write.csv(isangi_points,file.path(out,"Figure_1D_Isangi_positions.csv"),row.names=FALSE)
  saveRDS(list(A=pA,Bi=pBi,Bii=pBii,C=network_grob,D=pD,composite=figure),file.path(out,"Figure_1_editable_objects.rds"))
  caption <- paste0(
    "Figure 1. Resistance, genomic relatedness and environmental sampling of Salmonella in Blantyre, Malawi. ",
    "(A) Annual numbers of non-typhoidal Salmonella (NTS) blood-culture isolates at Queen Elizabeth Central Hospital ",
    "(QECH; bars, left axis), with the proportions in the recorded MDR-but-not-XDR and XDR categories (lines, right axis). ",
    "MDR denotes resistance to ampicillin, chloramphenicol and co-trimoxazole; XDR additionally includes ",
    "fluoroquinolone and third-generation cephalosporin resistance under the study's operational definition. ",
    "(B) SNP distances of 74 Malawian ST335 isolates from CAAP2A over the full collection period (i) and during ",
    "2020 (ii). Points are positioned at the midpoint of the recorded month; point area represents coincident ",
    "isolates with the same source and distance. CSF, cerebrospinal fluid. Colours in B and C identify isolate source. ",
    "(C) An undirected minimum spanning tree of the validated SNP matrix. Only identical complete distance profiles ",
    "are collapsed; zero-SNP edges are retained. Node area is proportional to isolate count; sectors show sources. ",
    "Bold edge labels give SNP distances; complete node membership is supplied in the accompanying source-data table. ",
    "Position and edge length are for display only and do not establish transmission direction. ",
    "(D) Environmental recovery locations of S. Isangi and other NTS, relative to QECH and the supplied waterway network. ",
    "Coordinates are not jittered; coincident observations overlap. The map is cropped to include all plotted sampling locations.")
  writeLines(caption,file.path(out,"Figure_1_caption.txt"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  input_paths <- na.omit(c(annual_path,time_path,matrix_path,bp,rp,ep,ip))
  utils::write.csv(data.frame(path=input_paths,md5=unname(tools::md5sum(input_paths))),
                   file.path(out,"input_checksums.csv"),row.names=FALSE)

  render <- function(path, dpi, tiff=TRUE) {
    temp <- paste0(path,".partial")
    if(tiff) ragg::agg_tiff(temp,width=FW,height=FH,units="mm",res=dpi,
                            compression="lzw",background="white") else
      ragg::agg_png(temp,width=FW,height=FH,units="mm",res=dpi,background="white")
    device <- grDevices::dev.cur()
    ok <- FALSE
    tryCatch({ print(figure); ok <- TRUE }, finally={ grDevices::dev.off(device) })
    if(!ok || !file.exists(temp) || file.info(temp)$size<=0) stop("Rendering failed: ",path)
    if(!file.rename(temp,path)) stop("Cannot finalise output: ",path)
  }
  tiff_path <- file.path(out,"Figure_1A_D_Genome_Medicine_600dpi.tiff")
  png_path <- file.path(out,"Figure_1A_D_preview.png")
  render(tiff_path,600,TRUE)
  render(png_path,180,FALSE)
  record("Saved TIFF: ",tiff_path)
  record("Canvas: 170 x 190 mm, 600 dpi. File size: ",round(file.info(tiff_path)$size/1024^2,2)," MB.")
  if(isTRUE(open_result) && unname(Sys.info()["sysname"])=="Darwin")
    try(system2("/usr/bin/open",c("-a","Preview",shQuote(tiff_path)),wait=FALSE),silent=TRUE)
  invisible(list(figure=figure,output_dir=out,tiff=tiff_path,preview=png_path))
}
