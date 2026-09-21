#!/usr/bin/env Rscript

ISANGI_VERSION <- "1.0.0"
ISANGI_HOME <- local({
  ff <- vapply(sys.frames(), function(e) if (!is.null(e$ofile)) as.character(e$ofile)[1L] else "", character(1))
  ff <- ff[nzchar(ff)]
  ca <- grep("^--file=", commandArgs(), value = TRUE)
  p <- if (length(ff)) tail(ff, 1L) else if (length(ca)) sub("^--file=", "", ca[1L]) else ""
  if (nzchar(p)) dirname(normalizePath(p, mustWork = FALSE)) else getwd()
})

isangi_config <- function(project_root = getwd(), input_root = file.path(ISANGI_HOME, "inputs")) {
  list(
    project_root = path.expand(project_root), input_root = path.expand(input_root),
    source_zip = character(), search_roots = character(), paths = list(),

    output_parent = "outputs/publication_figures", output_dir = NULL,
    targets = "all", formats = c("tiff", "pdf", "svg", "png"),
    dpi = 600L, preview_dpi = 160L, seed = 20260919L,
    font_family = "sans", font_scale = 1, branch_units = "model branch-length units",
    fonts = list(legend = 7.6, legend_heading = 8.0, tip = 7.6,
                 heat_column = 7.2, matrix_label = 7.0, axis = 8.4, cell = 7.4),

    sizes = list(),
    save_grobs = TRUE, return_figures = TRUE, export_long_matrices = FALSE,
    create_assemblies = TRUE, fail_at_end = TRUE,

    check_frozen_results = TRUE,
    expected = list(master = 345L, cc25 = 327L, st335 = 224L, phenotype = 42L,
                    isangi_phenotype = 19L, comparator = 23L, ena = 80L),
    palette = list(
      country = c(Malawi = "#CD414B", `South Africa` = "#287DA8", Nigeria = "#53A058",
        `United Kingdom` = "#87559D", Brazil = "#E89032", `United States` = "#D8BC47",
        Mexico = "#91663D", Taiwan = "#D979AF", Mozambique = "#49B7B5", Uganda = "#A75F87", Other = "#BFC3C9"),
      source = c(`River water` = "#69B99C", `QECH clinical` = "#EE996D",
        `Neonatal unit sampling` = "#9CACD0", `South African outbreak` = "#D58CAF", Other = "#F1F2F4"),
      timeline = c(`Neonatal unit sampling` = "#9CACD0", `River water` = "#69B99C",
                   `Blood culture` = "#EE996D", CSF = "#74528D"),
      year = c(`2000–2004` = "#F0F7EB", `2005–2008` = "#CDE6C9", `2009–2012` = "#A4D8C5",
        `2013–2016` = "#76C3C6", `2017–2020` = "#43A4C4", `2021–2024` = "#206B9A", Other = "#D6D8DD"),
      lineage = c(`1` = "#F2DBA9", `2` = "#CEE5C6", `3` = "#DFCCE7", `4` = "#C3D8E8", `5` = "#EDC4C8", Other = "#E5E5E5"),
      st = c(`216` = "#E7D97C", `335` = "#62B2AB", `5028` = "#A790BC", Other = "#BBBBBB"),
      amr_count = c(`0` = "#F2F0F7", `1–3` = "#CECAE1", `4–6` = "#A7A0C8", `7–10` = "#7966AB", `11–13` = "#4D2C83", Other = "#BFC3C9"),
      phenotype = c(S = "#DAECE6", I = "#E8BC56", R = "#C94E55", Other = "#F1F2F4"),
      bsi = c(XDR = "#C94E55", MDR = "#E7A355", Other = "#B4CCD8"),
      gene = c(`Not detected` = "#E2E5E9", Detected = "#599EAF", Other = "#F1F2F4"),
      status = c(`Other profile` = "#E2E5E9", `Profile / marker` = "#C94E55", Other = "#F1F2F4"),
      ink = "#30373F", grid = "#DCE1E5",

      snp = c("#FDE725", "#B5DE2B", "#6DCD59", "#35B779", "#1F9E89",
              "#26828E", "#31688E", "#3E4989", "#482878", "#440154")
    ),
    circular = list(gap_degrees = 18, clockwise = TRUE, tip_diameter_mm = 1.30,
      branch_lwd = 0.55, centre_fraction = c(0.372, 0.505), radius_fraction = 0.452,
      inner_radius = 0.08, tree_radius = 0.78, shade_radius = 0.812,
      ring_radii = c(0.866, 0.940, 1.014, 1.088), ring_width = 0.050,
      ring_labels = c("Country", "Source", "Year", "AMR classes"),
      ring_label_font = 7.4, scale_length = 0.02),
    tree = list(main_tip_diameter_mm = 1.2, detail_tip_diameter_mm = 1.5,
      matched_tip_diameter_mm = 2.0, branch_lwd = 0.50,

      main_tree_box = c(.025, .185, .250, .770),
      main_heat_box = c(.295, .185, .515, .770),
      main_legend_box = c(.840, .12, .150, .835),
      main_tracks = "Year", detail_cc25_tracks = c("Country", "Source", "Year", "Lineage", "AMR classes", "XDR profile")),
    genes = list(main = NULL, supplementary = NULL, genes_per_block = 28L),
    panels = list(use_frozen_membership = TRUE, max_rows = 34L),
    phenotype = list(cell_text = "zone",

      resistant_codes = "R", mdr_agents = c("AMP10", "C30", "SXT25"),
      fluoroquinolone_agent = "PEF5", cephalosporin_agents = c("CPD10", "CTX5")),
    matrix = list(transform = "log1p",
      max_cc25 = 12000, max_st335 = 560, draw_cells_as_vectors = TRUE,
      metadata = c("Country", "Year", "Source", "Lineage", "XDR profile"),
      detail_labels = "isolate_accession", overview_tick_every = 50L),
    timeline = list(reference = "CAAP2A",
      point_min_mm = 1.8, point_count_mm2 = 0.40, first_year = 2018, last_year = 2023),
    network = list(iterations = 450L, collision_iterations = 300L,
      k = 0.25, base_node_radius = .025, label_edges = TRUE, label_node_counts = TRUE,
      coordinates = NULL, collapse = "zero_components"),

    legacy = list(mode = "retain",
      map_gpkg = NULL, map_membership = NULL, map_epsg = 32736L,
      ncbi_counts = NULL, plasmids = NULL, plasmid_features = NULL,
      plasmid_regions = NULL, plasmid_alignments = NULL,
      plasmid_nodes = NULL, plasmid_edges = NULL, plasmid_positions = NULL,
      confocal_manifest = NULL, biofilm_values = NULL, biofilm_annotations = NULL,
      mouse_survival = NULL, allow_title_masks = TRUE),

    custom_builders = list()
  )
}

ISANGI_AGENTS <- c("AMP10", "AZM15", "CIP5", "FOX30", "MEM10", "C30", "PEF5",
  "FEP30", "SXT25", "TCG15", "IPM10", "AK30", "CPD10", "ATM30", "CN10", "TZP36", "CTX5")
ISANGI_EXTRAS <- c("ERR9939684", "ERR11201713", "ERR9939690", "SRR3049609")

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
other <- function(x) {
  x <- trimws(as.character(x)); x[blank(x) | tolower(x) %in% c("unknown", "not known", "unavailable", "na", "n/a")] <- "Other"; x
}
clean_id <- function(x) {
  x <- trimws(as.character(x)); x <- sub("_assembly_contigs\\.fasta$", "", x)
  x <- sub("_S[0-9]+_L[0-9]+$", "", x); gsub("#", "_", x, fixed = TRUE)
}
unique_ids <- function(x, what) {
  assert(!any(blank(x)), paste(what, "has blank IDs"))
  assert(!anyDuplicated(x), paste(what, "has duplicate or normalisation-colliding IDs"))
}
need_cols <- function(x, cols, what) {
  assert(all(cols %in% names(x)), paste(what, "is missing:", paste(setdiff(cols, names(x)), collapse = ", ")))
}
number <- function(x, what, allow_missing = FALSE) {
  z <- suppressWarnings(as.numeric(x))
  assert(!any(!blank(x) & is.na(z)), paste("Non-numeric value in", what))
  if (!allow_missing) assert(all(is.finite(z)), paste("Missing/non-finite value in", what))
  z
}
read_tab <- function(path) {
  sep <- if (grepl("\\.(tsv|tab)$", path, ignore.case = TRUE)) "\t" else ","
  z <- utils::read.table(path, header = TRUE, sep = sep, quote = '"', comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = c("", "NA"), fill = FALSE)
  assert(!anyDuplicated(names(z)), paste("Duplicate columns:", path)); z
}
write_tab <- function(x, file) utils::write.csv(x, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
canonical_filename <- function(x) {
  x <- basename(x)
  x <- gsub("(\\s*\\([0-9]+\\))+(?=\\.[^.]+$|$)", "", x, perl = TRUE)
  tolower(trimws(x))
}
index_files <- function(roots) {
  roots <- unique(path.expand(roots[!blank(roots)])); roots <- roots[dir.exists(roots)]
  ff <- unique(unlist(lapply(roots, list.files, recursive = TRUE, full.names = TRUE, all.files = FALSE)))
  ff <- ff[!dir.exists(ff) & !grepl("/__MACOSX/|/\\.Rproj\\.user/|/\\.git/", ff)]
  data.frame(path = ff, basename = canonical_filename(ff), stringsAsFactors = FALSE)
}
extract_sources <- function(zips) {
  out <- character()
  for (z in zips) {
    z <- path.expand(z); assert(file.exists(z), paste("Source ZIP not found:", z))
    a <- utils::unzip(z, list = TRUE)$Name; a <- gsub("\\\\", "/", a)
    assert(!any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", a)), "Unsafe path in source ZIP")
    dest <- tempfile("isangi_sources_"); dir.create(dest)
    keep <- a[!grepl("(^|/)__MACOSX/|(^|/)\\._", a)]
    utils::unzip(z, files = keep, exdir = dest); out <- c(out, dest)
  }; out
}
resolve_input <- function(ctx, role, basenames, required = TRUE) {
  override <- ctx$cfg$paths[[role]]
  if (!is.null(override)) {
    p <- path.expand(override); assert(file.exists(p), paste("Configured input absent:", role, p))
  } else {
    p <- NULL

    for (idx in list(ctx$bundled, ctx$external)) {
      for (b in basenames) {
        q <- idx$path[idx$basename == canonical_filename(b)]
        if (length(q)) {
          hashes <- unname(tools::md5sum(q))
          assert(length(unique(hashes)) == 1L, paste("Conflicting copies for", role,
            "— set cfg$paths explicitly:\n", paste(q, collapse = "\n")))
          p <- sort(q)[1L]; break
        }
      }
      if (!is.null(p)) break
    }
  }
  if (is.null(p)) {
    if (required) stop("Missing ", role, ": ", paste(basenames, collapse = " / "), call. = FALSE)
    return(NULL)
  }
  ctx$used[[role]] <- normalizePath(p, mustWork = TRUE); p
}
validate_matrix <- function(D, label = "SNP matrix") {
  assert(is.matrix(D) && is.numeric(D) && nrow(D) == ncol(D), paste(label, "is not square/numeric"))
  rownames(D) <- clean_id(rownames(D)); colnames(D) <- clean_id(colnames(D))
  unique_ids(rownames(D), paste(label, "rows")); unique_ids(colnames(D), paste(label, "columns"))
  assert(setequal(rownames(D), colnames(D)), paste(label, "row and column sets differ"))
  D <- D[, match(rownames(D), colnames(D)), drop = FALSE]
  assert(all(is.finite(D)) && all(D >= 0) && all(D == round(D)), paste(label, "has missing/negative/non-integer counts"))
  assert(all(diag(D) == 0) && all(D == t(D)), paste(label, "is asymmetric or has a nonzero diagonal")); D
}
read_matrix <- function(p) {
  z <- read_tab(p); assert(ncol(z) == nrow(z) + 1L, paste("Expected one ID field plus a square matrix:", p))
  ids <- clean_id(z[[1L]]); m <- do.call(cbind, lapply(z[-1L], number, what = p))
  dimnames(m) <- list(ids, clean_id(names(z)[-1L])); validate_matrix(m, basename(p))
}
export_matrix <- function(D, name, out, long = FALSE) {
  write_tab(data.frame(isolate_id = rownames(D), D, check.names = FALSE), file.path(out, paste0(name, ".csv")))
  if (long) {
    ij <- which(upper.tri(D), arr.ind = TRUE)
    write_tab(data.frame(isolate_a = rownames(D)[ij[, 1]], isolate_b = colnames(D)[ij[, 2]],
      snp_distance = D[ij]), file.path(out, paste0(name, "_pairs.csv")))
  }
}
check_dependencies <- function(formats = c("tiff", "png", "svg")) {
  pkgs <- c("ape", "jsonlite", "png")
  if (any(c("tiff", "png") %in% formats)) pkgs <- c(pkgs, "ragg")
  if ("svg" %in% formats) pkgs <- c(pkgs, "svglite")
  absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  assert(!length(absent), paste("Missing R packages:", paste(absent, collapse = ", "),
    "\nRun install_isangi_dependencies.R once. Rendering never installs or updates packages."))
}

tree_structure <- function(tr) {
  nt <- length(tr$tip.label); nn <- nt + tr$Nnode
  tr$tip.label <- clean_id(tr$tip.label); unique_ids(tr$tip.label, "Tree")
  assert(!is.null(tr$edge.length) && length(tr$edge.length) == nrow(tr$edge), "Tree lacks branch lengths")
  assert(all(is.finite(tr$edge.length)) && all(tr$edge.length >= 0), "Tree has invalid branch lengths")
  root <- unique(setdiff(tr$edge[, 1], tr$edge[, 2])); assert(length(root) == 1L, "Tree root is not unique")
  children <- split(tr$edge[, 2], tr$edge[, 1]); edge_len <- numeric(nn)
  edge_len[tr$edge[, 2]] <- tr$edge.length
  depth <- numeric(nn); order <- integer()
  visit <- function(node) {
    if (node <= nt) { order <<- c(order, node); return(invisible(NULL)) }
    for (ch in children[[as.character(node)]]) { depth[ch] <<- depth[node] + edge_len[ch]; visit(ch) }
  }; visit(root)
  list(tree = tr, nt = nt, nn = nn, root = root, children = children,
       depth = depth, order = order, ids = tr$tip.label[order])
}
layout_tree <- function(st, ids = st$ids) {
  assert(all(ids %in% st$ids) && !anyDuplicated(ids) && length(ids) > 0L, "Invalid tree subset")
  ids <- st$ids[st$ids %in% ids]; active <- logical(st$nn)
  descend <- function(v) {
    active[v] <<- if (v <= st$nt) st$tree$tip.label[v] %in% ids else any(vapply(st$children[[as.character(v)]], descend, logical(1)))
    active[v]
  }; descend(st$root)
  root <- st$root
  while (root > st$nt) {
    ch <- st$children[[as.character(root)]]; ch <- ch[active[ch]]
    if (length(ch) != 1L) break
    root <- ch[1L]
  }
  y <- rep(NA_real_, st$nn); ix <- match(ids, st$tree$tip.label); y[ix] <- seq_along(ids)
  yy <- function(v) {
    if (!is.na(y[v])) return(y[v])
    ch <- st$children[[as.character(v)]]; ch <- ch[active[ch]]
    z <- vapply(ch, yy, numeric(1)); y[v] <<- (z[1L] + tail(z, 1L)) / 2; y[v]
  }; yy(root)
  x <- st$depth - st$depth[root]; seg <- list(); k <- 0L
  traverse <- function(v) {
    if (v <= st$nt) return(invisible(NULL))
    ch <- st$children[[as.character(v)]]; ch <- ch[active[ch]]
    k <<- k + 1L; seg[[k]] <<- c(x[v], y[ch[1L]], x[v], y[tail(ch, 1L)])
    for (u in ch) { k <<- k + 1L; seg[[k]] <<- c(x[v], y[u], x[u], y[u]); traverse(u) }
  }; traverse(root)
  sm <- if (length(seg)) do.call(rbind, seg) else matrix(numeric(), ncol = 4L)
  colnames(sm) <- c("x0", "y0", "x1", "y1")
  list(ids = ids, root = root, x = x, y = y, segments = sm, nodes = ix,
       tip_x = x[ix], tip_y = seq_along(ids), max_x = max(c(x[ix], 1e-12)), structure = st)
}

classify_phenotypes <- function(p, cfg) {
  R <- as.matrix(p[, ISANGI_AGENTS, drop = FALSE]) %in% cfg$phenotype$resistant_codes
  dim(R) <- c(nrow(p), length(ISANGI_AGENTS)); dimnames(R) <- list(p$isolate_id, ISANGI_AGENTS)

  R[is.na(as.matrix(p[, ISANGI_AGENTS, drop = FALSE]))] <- NA
  all_row <- function(m) apply(m, 1L, all)
  any_row <- function(m) apply(m, 1L, any)
  mdr <- all_row(R[, cfg$phenotype$mdr_agents, drop = FALSE])
  xdr <- mdr & R[, cfg$phenotype$fluoroquinolone_agent] & any_row(R[, cfg$phenotype$cephalosporin_agents, drop = FALSE])
  data.frame(isolate_id = p$isolate_id, group = p$group, MDR_phenotype = as.integer(mdr),
             XDR_phenotype = as.integer(xdr), recorded_R_count = rowSums(R), stringsAsFactors = FALSE)
}
classify_genomes <- function(g, markers) {
  out <- data.frame(isolate_id = rownames(g), check.names = FALSE)
  for (cat in names(markers)) {
    gg <- unlist(markers[[cat]], use.names = FALSE)
    assert(all(gg %in% colnames(g)), paste("Missing marker columns for", cat))
    out[[cat]] <- as.integer(rowSums(g[, gg, drop = FALSE]) > 0)
  }
  mdr <- out$Penicillin == 1 & out$Chloramphenicol == 1 & out$Trimethoprim == 1 & out$Sulfonamide == 1

  out$MDR_genomic_profile <- as.integer(mdr)
  out$XDR_genomic_profile <- as.integer(mdr & out$Fluoroquinolone == 1 & out[["Third-generation cephalosporin"]] == 1)
  out$CARB_R_genomic_marker <- out$Carbapenem; out$MAC_R_genomic_marker <- out$Macrolide
  out
}
make_panels <- function(st, meta, population, max_rows) {
  line <- other(meta[st$ids, "lineage"]); rr <- rle(line); start <- cumsum(c(1L, head(rr$lengths, -1L)))
  ix <- list(); jj <- 0L
  for (i in seq_along(rr$lengths)) {
    pos <- seq.int(start[i], length.out = rr$lengths[i])

    groups <- split(pos, ceiling(seq_along(pos) / max_rows))
    for (pp in groups) {
      jj <- jj + 1L; pid <- sprintf("%s-%02d", population, jj)
      ix[[jj]] <- data.frame(population = population, panel_id = pid, isolate_id = st$ids[pp],
        population_tip_order = pp, panel_tip_order = seq_along(pp), lineage = line[pp],
        run_accession = meta[st$ids[pp], "run_accession"], stringsAsFactors = FALSE)
    }
  }; do.call(rbind, ix)
}
panel_summary <- function(pi) {
  do.call(rbind, lapply(unique(pi$panel_id), function(pid) {
    p <- pi[pi$panel_id == pid, , drop = FALSE]
    data.frame(population = p$population[1], panel_id = pid, lineage = p$lineage[1], n_isolates = nrow(p),
      first_tip_position = min(as.integer(p$population_tip_order)), last_tip_position = max(as.integer(p$population_tip_order)),
      first_isolate = p$isolate_id[1], last_isolate = tail(p$isolate_id, 1),
      partition_basis = "Consecutive tip-order subset within recorded lineage", stringsAsFactors = FALSE)
  }))
}
load_isangi <- function(cfg, output_dir = NULL) {
  ctx <- new.env(parent = emptyenv()); ctx$cfg <- cfg; ctx$used <- list()
  extra_roots <- extract_sources(cfg$source_zip)
  ctx$bundled <- index_files(cfg$input_root)
  ctx$external <- index_files(c(cfg$search_roots, extra_roots, cfg$project_root))
  tab <- function(role, files, required = TRUE) {
    p <- resolve_input(ctx, role, files, required); if (is.null(p)) NULL else read_tab(p)
  }
  sc <- jsonlite::fromJSON(resolve_input(ctx, "scientific_config", c("scientific_config.json", "config.json")), simplifyVector = TRUE)
  assert(identical(as.character(sc$antibiotics), ISANGI_AGENTS), "Scientific config does not contain the exact 17 recorded agents")
  m <- tab("metadata", "isolate_accession_line_list_345.csv"); need_cols(m, c("isolate_id", "run_accession", "sample_accession", "country", "year", "lineage", "sequence_type", "source_group", "AMRFinder_class_count"), "Master metadata")
  m$isolate_id <- clean_id(m$isolate_id); unique_ids(m$isolate_id, "Master metadata"); rownames(m) <- m$isolate_id
  if (!"plot_label" %in% names(m)) m$plot_label <- ifelse(m$isolate_id == m$run_accession, m$isolate_id, paste(m$isolate_id, m$run_accession, sep = " | "))

  ct <- tree_structure(ape::read.tree(resolve_input(ctx, "cc25_tree", c("CC25_tree.nwk", "Supplement_CC25_plotted.nwk"))))
  st <- tree_structure(ape::read.tree(resolve_input(ctx, "st335_tree", c("ST335_tree.nwk", "Supplement_ST335_plotted.nwk"))))
  assert(all(ct$ids %in% m$isolate_id) && all(st$ids %in% ct$ids), "Tree IDs do not reconcile to the master cohort")
  full_path <- resolve_input(ctx, "complete_matrix", c("CC25_complete_331.csv", "ebg_25_snpmatrix.tsv"), FALSE)
  full <- if (!is.null(full_path)) read_matrix(full_path) else NULL
  cp <- resolve_input(ctx, "aligned_matrix", c("CC25_tree_aligned_327.csv", "CC25_pairwise_SNP_matrix_tree_order.csv"), is.null(full))
  old_cc <- if (!is.null(cp)) read_matrix(cp) else NULL
  source_matrix <- full %||% old_cc
  assert(all(ct$ids %in% rownames(source_matrix)), "CC25 tree tips are absent from the supplied numeric matrix")
  cc <- source_matrix[ct$ids, ct$ids, drop = FALSE]; ss <- cc[st$ids, st$ids, drop = FALSE]
  if (!is.null(old_cc)) assert(all(old_cc[ct$ids, ct$ids] == cc), "Complete and previously aligned CC25 values conflict")
  sp <- resolve_input(ctx, "st335_matrix", c("ST335_tree_aligned_224.csv", "ST335_pairwise_SNP_matrix_tree_order.csv"), FALSE)
  if (!is.null(sp)) { old_st <- read_matrix(sp); assert(all(old_st[st$ids, st$ids] == ss), "ST335 matrix conflicts with its CC25 subset") }
  if (!is.null(full)) {
    extra <- setdiff(rownames(full), ct$ids)
    if (cfg$check_frozen_results) assert(setequal(extra, ISANGI_EXTRAS), "Unexpected matrix-only isolate set")
    extra <- sort(extra); full <- full[c(ct$ids, extra), c(ct$ids, extra), drop = FALSE]
  }
  reconcile <- data.frame(isolate_id = ISANGI_EXTRAS, in_CC25_tree = FALSE,
    complete_pairwise_distances_supplied = !is.null(full),
    treatment = if (!is.null(full)) "Retained in complete matrix; not added to tree or study cohort" else
      "Matrix-source isolate; complete pairwise distances not supplied", stringsAsFactors = FALSE)
  p <- tab("phenotype_calls", "phenotype_categories_as_recorded_42.csv")
  z <- tab("phenotype_zones", "phenotype_zone_diameters_42.csv")
  need_cols(p, c("isolate_id", "group", ISANGI_AGENTS), "Recorded AST"); need_cols(z, c("isolate_id", "group", ISANGI_AGENTS), "Recorded zones")
  p$isolate_id <- clean_id(p$isolate_id); z$isolate_id <- clean_id(z$isolate_id)
  unique_ids(p$isolate_id, "AST"); unique_ids(z$isolate_id, "Zone data")
  assert(setequal(p$isolate_id, z$isolate_id), "AST/zone isolate sets differ")
  z <- z[match(p$isolate_id, z$isolate_id), , drop = FALSE]
  assert(identical(p$group, z$group), "AST/zone group labels conflict")
  assert(all(as.matrix(p[, ISANGI_AGENTS]) %in% c("S", "I", "R")), "Unrecognised or missing recorded AST category")
  for (a in ISANGI_AGENTS) { z[[a]] <- number(z[[a]], a); assert(all(z[[a]] >= 0), "Negative zone diameter") }
  rownames(p) <- p$isolate_id; rownames(z) <- z$isolate_id
  assert(all(p$isolate_id[p$group == "Isangi"] %in% st$ids), "Tested Isangi genomes lack an exact ST335 match")
  ph <- classify_phenotypes(p, cfg)
  gd <- tab("determinants", "genomic_determinants_345.csv"); unique_ids(clean_id(gd$isolate_id), "Genomic determinants")
  genes <- as.character(sc$gene_order); need_cols(gd, c("isolate_id", genes), "Genomic calls")
  assert(setequal(clean_id(gd$isolate_id), m$isolate_id), "Genomic calls and master isolate sets differ")
  gd <- gd[match(m$isolate_id, clean_id(gd$isolate_id)), , drop = FALSE]
  g <- do.call(cbind, lapply(gd[, genes, drop = FALSE], number, what = "Genomic call"))
  dimnames(g) <- list(m$isolate_id, genes); assert(all(g %in% c(0, 1)), "Genomic calls must be explicit zero or one")
  gf <- classify_genomes(g, sc$genomic_category_markers); rownames(gf) <- gf$isolate_id
  for (col in c("MDR_genomic_profile", "XDR_genomic_profile", "CARB_R_genomic_marker", "MAC_R_genomic_marker")) {
    if (cfg$check_frozen_results) assert(all(number(m[[col]], col) == gf[m$isolate_id, col]), paste("Genomic rule differs from frozen", col))
    m[[col]] <- gf[m$isolate_id, col]
  }
  main_genes <- cfg$genes$main %||% as.character(sc$main_gene_order)
  supp_genes <- cfg$genes$supplementary %||% genes
  assert(!anyDuplicated(main_genes) && all(main_genes %in% genes), "Invalid main gene order")
  assert(setequal(supp_genes, genes) && !anyDuplicated(supp_genes), "Supplementary genes must retain all curated determinant columns")
  assert(all(colnames(g)[colSums(g[st$ids, , drop = FALSE]) > 0] %in% main_genes), "Main gene list drops a determinant observed in ST335")
  pi <- if (cfg$panels$use_frozen_membership) tab("panels", "panel_isolate_index.csv") else
    rbind(make_panels(ct, m, "CC25", cfg$panels$max_rows), make_panels(st, m, "ST335", cfg$panels$max_rows))
  for (pop in c("CC25", "ST335")) {
    ids <- if (pop == "CC25") ct$ids else st$ids; pp <- pi[pi$population == pop, , drop = FALSE]
    unique_ids(pp$isolate_id, paste(pop, "panel index")); assert(setequal(ids, pp$isolate_id), paste(pop, "panel coverage mismatch"))
    for (pid in unique(pp$panel_id)) {
      q <- pp[pp$panel_id == pid, , drop = FALSE]; oi <- match(q$isolate_id, ids)
      assert(all(diff(oi) == 1) && length(unique(other(m[q$isolate_id, "lineage"]))) == 1L, paste("Non-consecutive or mixed-lineage panel", pid))
    }
  }
  ps <- panel_summary(pi)
  annual <- tab("bsi_annual", "BSI_annual_resistance_summary_2015_2023.csv")
  need_cols(annual, c("year", "resistance_group", "n", "percent"), "BSI summary")
  annual$year <- number(annual$year, "BSI year"); annual$n <- number(annual$n, "BSI count")
  annual$percent <- number(annual$percent, "BSI percent")
  assert(all(annual$resistance_group %in% c("XDR", "MDR", "Other")), "Unexpected BSI group")
  for (yr in unique(annual$year)) {
    q <- annual[annual$year == yr, ]; assert(abs(sum(q$percent) - 100) < 1e-7, "BSI percentages do not sum to 100")
    assert(all(abs(q$percent - 100 * q$n / sum(q$n)) < 1e-7), "BSI counts and percentages disagree")
  }
  timeline <- tab("timeline", "Malawi_timeline_74_month_resolution.csv")
  need_cols(timeline, c("isolate_id", "collection_month", "source", "SNP_distance_to_reference"), "Timeline")
  timeline$isolate_id <- clean_id(timeline$isolate_id); unique_ids(timeline$isolate_id, "Timeline")
  assert(all(timeline$isolate_id %in% st$ids) && cfg$timeline$reference %in% st$ids, "Timeline/ref IDs do not match ST335")
  assert(all(grepl("^[0-9]{4}-[0-9]{2}$", timeline$collection_month)), "Timeline requires explicit month-resolution dates")

  timeline$plot_date <- as.Date(paste0(timeline$collection_month, "-15"))
  assert(!anyNA(timeline$plot_date), "Invalid collection month")
  timeline$SNP_distance_to_reference <- as.numeric(ss[timeline$isolate_id, cfg$timeline$reference])
  timeline$reference_isolate <- cfg$timeline$reference
  mw <- ss[timeline$isolate_id, timeline$isolate_id, drop = FALSE]

  mp <- resolve_input(ctx, "malawi_matrix", "Malawi_74_source_matrix.csv", FALSE)
  if (!is.null(mp)) {
    old <- read_matrix(mp); assert(setequal(rownames(old), timeline$isolate_id), "Malawi matrix membership conflict")
    assert(all(old == ss[rownames(old), colnames(old)]), "Malawi source matrix differs from validated ST335 values")
    mw <- old
  }
  extras_meta <- tab("matrix_metadata", "CC25_331_matrix_isolate_index.csv", FALSE)
  if (!is.null(extras_meta)) { extras_meta$isolate_id <- clean_id(extras_meta$isolate_id); unique_ids(extras_meta$isolate_id, "Full-matrix metadata"); rownames(extras_meta) <- extras_meta$isolate_id }
  dict <- tab("antibiotic_dictionary", "antibiotic_dictionary_17.csv")
  assert(identical(as.character(dict$code), ISANGI_AGENTS), "Antibiotic dictionary differs from recorded columns")
  summ <- do.call(rbind, lapply(ISANGI_AGENTS, function(a) {
    x <- p[p$group == "Isangi", a]; y <- p[p$group == "Comparator", a]
    r1 <- sum(x == "R"); r2 <- sum(y == "R")
    data.frame(code = a, antibiotic = dict$antibiotic[match(a, dict$code)],
      Isangi_S = sum(x == "S"), Isangi_I = sum(x == "I"), Isangi_R = r1, Isangi_n = length(x), Isangi_R_percent = 100 * r1 / length(x),
      Comparator_S = sum(y == "S"), Comparator_I = sum(y == "I"), Comparator_R = r2, Comparator_n = length(y), Comparator_R_percent = 100 * r2 / length(y),
      Fisher_exact_p = stats::fisher.test(matrix(c(r1, length(x)-r1, r2, length(y)-r2), 2L, byrow = TRUE))$p.value,
      stringsAsFactors = FALSE)
  }))
  ena <- tab("ena", "ENA_deposited_isolates_80.csv")
  assert(all(ena$isolate_id %in% m$isolate_id) && all(m[ena$isolate_id, "run_accession"] == ena$run_accession), "Deposited accession crosswalk mismatch")
  checks <- c(master_unique = TRUE, matrix_valid = TRUE, cc25_st335_exact_subset = TRUE,
    exact_17_agents = TRUE, measured_phenotype_not_reinterpreted = TRUE, panel_coverage = TRUE,
    explicit_matrix_only_reconciliation = TRUE, ena_exact_join = TRUE)
  if (cfg$check_frozen_results) {
    assert(nrow(m) == cfg$expected$master && length(ct$ids) == cfg$expected$cc25 && length(st$ids) == cfg$expected$st335, "Frozen genome counts changed")
    assert(nrow(p) == 42L && sum(p$group == "Isangi") == 19L && sum(p$group == "Comparator") == 23L, "Frozen phenotype group counts changed")
    assert(sum(ph$XDR_phenotype[ph$group == "Isangi"]) == 12L, "Frozen phenotypic XDR count changed")
    assert(sum(gf[st$ids, "XDR_genomic_profile"]) == 199L, "Frozen ST335 genomic XDR count changed")
    assert(sum(gf[st$ids, "CARB_R_genomic_marker"]) == 6L, "Frozen carbapenemase marker count changed")
    assert(all(gf[p$isolate_id[p$group == "Isangi"], "XDR_genomic_profile"] == 1L), "Matched genome profiles changed")
    assert(nrow(ena) == cfg$expected$ena && all(!blank(m$run_accession)) && all(!blank(m$sample_accession)), "Frozen accession coverage changed")
    old_summary <- tab("phenotype_summary", "phenotype_summary_17.csv")
    old_summary <- old_summary[match(summ$code, old_summary$code), ]
    assert(all(summ$Isangi_R == number(old_summary$Isangi_R, "Isangi R")) &&
      all(summ$Comparator_R == number(old_summary$Comparator_R, "Comparator R")) &&
      max(abs(summ$Fisher_exact_p - number(old_summary$Fisher_exact_p, "Fisher P"))) < 1e-10, "AST summary differs from frozen results")
    checks <- c(checks, frozen_345_327_224 = TRUE, frozen_12_of_19 = TRUE, frozen_199_of_224 = TRUE, frozen_6_carbapenemase = TRUE)
  }
  legacy <- jsonlite::fromJSON(resolve_input(ctx, "legacy_manifest", "legacy_figure_manifest.json"), simplifyVector = FALSE)
  for (j in seq_along(legacy)) legacy[[j]]$path <- resolve_input(ctx, paste0("legacy_", j), legacy[[j]]$file)
  cap <- list()
  for (role in c("main_figure_captions.json", "supplementary_figure_captions.json")) {
    fp <- resolve_input(ctx, role, role, FALSE); if (!is.null(fp)) cap <- c(cap, jsonlite::fromJSON(fp, simplifyVector = FALSE))
  }
  registry <- jsonlite::fromJSON(resolve_input(ctx, "figure_registry", "figure_registry.json"), simplifyVector = TRUE)
  dat <- list(meta = m, genes = g, genomic_flags = gf, scientific_config = sc,
    cc_tree = ct, st_tree = st, cc = cc, st = ss, full = full, malawi = mw,
    phenotype = p, zones = z, phenotype_flags = ph, phenotype_summary = summ,
    antibiotic_dictionary = dict, annual = annual, timeline = timeline,
    panels = pi, panel_summary = ps, extras_meta = extras_meta, reconciliation = reconcile,
    main_genes = main_genes, supp_genes = supp_genes, legacy = legacy,
    captions = cap, registry = registry, input_paths = ctx$used, checks = checks)
  if (!is.null(output_dir)) export_analysis(dat, cfg, output_dir)
  dat
}
export_analysis <- function(d, cfg, out) {
  for (sub in c("source_data", "matrices", "audit", "captions", "objects", "main", "supplementary", "assemblies", "previews", "vector")) dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
  sd <- file.path(out, "source_data")
  for (pair in list(c("meta", "isolate_accession_line_list_345"), c("genomic_flags", "genomic_category_flags_R"),
    c("phenotype", "phenotype_categories_as_recorded_42"), c("zones", "phenotype_zone_diameters_42"),
    c("phenotype_flags", "phenotype_isolate_flags_R"), c("phenotype_summary", "phenotype_summary_17_R"),
    c("antibiotic_dictionary", "antibiotic_dictionary_17"), c("panels", "panel_isolate_index"),
    c("panel_summary", "panel_summary"), c("reconciliation", "CC25_matrix_tree_reconciliation_4"),
    c("annual", "BSI_annual_resistance_summary_2015_2023"))) write_tab(d[[pair[1]]], file.path(sd, paste0(pair[2], ".csv")))
  tt <- d$timeline; tt$plot_date <- NULL; write_tab(tt, file.path(sd, "Malawi_timeline_month_resolution_R.csv"))
  write_tab(data.frame(isolate_id = rownames(d$genes), d$genes, check.names = FALSE), file.path(sd, "genomic_determinants_345.csv"))
  for (p in list(c("cc", "CC25_tree_aligned_327"), c("st", "ST335_tree_aligned_224"), c("malawi", "Malawi_74_source_matrix"))) export_matrix(d[[p[1]]], p[2], file.path(out, "matrices"), cfg$export_long_matrices)
  if (!is.null(d$full)) export_matrix(d$full, "CC25_complete_331", file.path(out, "matrices"), cfg$export_long_matrices)
  for (pid in unique(d$panels$panel_id)) {
    pp <- d$panels[d$panels$panel_id == pid, ]; mm <- if (pp$population[1] == "CC25") d$cc else d$st
    export_matrix(mm[pp$isolate_id, pp$isolate_id, drop = FALSE], pid, file.path(out, "matrices"), cfg$export_long_matrices)
  }
  for (k in c("cc", "st", "full")) if (!is.null(d[[k]])) {
    ids <- rownames(d[[k]]); q <- metadata_for(d, ids)
    write_tab(data.frame(matrix_position = seq_along(ids), q, check.names = FALSE), file.path(sd, paste0(k, "_matrix_index_R.csv")))
  }
  pp <- unlist(d$input_paths, use.names = TRUE)
  write_tab(data.frame(role = names(pp), resolved_path = unname(pp), md5 = unname(tools::md5sum(pp))), file.path(out, "audit", "resolved_inputs.csv"))
  write_tab(data.frame(check = names(d$checks), result = ifelse(d$checks, "PASS", "FAIL")), file.path(out, "audit", "input_checks.csv"))
}

native_page <- function(w, h, cfg) {
  grid::pushViewport(grid::viewport(xscale = c(0, w), yscale = c(0, h), clip = "off",
    gp = grid::gpar(fontfamily = cfg$font_family)))
}
text_mm <- function(label, x, y, cfg, size = cfg$fonts$axis, colour = cfg$palette$ink,
                    just = "left", rot = 0, bold = FALSE, name = NULL) {
  grid::grid.text(label, x = x, y = y, default.units = "native", just = just, rot = rot,
    name = name, gp = grid::gpar(fontsize = size * cfg$font_scale, col = colour,
      fontfamily = cfg$font_family, fontface = if (bold) "bold" else "plain", lineheight = 1.15))
}
rect_mm <- function(x, y, width, height, fill, border = NA, lwd = .4, name = NULL) {
  grid::grid.rect(x = x, y = y, width = width, height = height, default.units = "native", name = name,
    gp = grid::gpar(fill = fill, col = border, lwd = lwd))
}
line_mm <- function(x0, y0, x1, y1, colour = "#30373F", lwd = .6, name = NULL) {
  grid::grid.segments(x0, y0, x1, y1, default.units = "native", name = name,
    gp = grid::gpar(col = colour, lwd = lwd, lineend = "butt"))
}
poly_mm <- function(x, y, fill = NA, border = NA, lwd = .4, name = NULL) {
  grid::grid.polygon(x, y, default.units = "native", name = name,
    gp = grid::gpar(fill = fill, col = border, lwd = lwd))
}
point_mm <- function(x, y, fill, diameter = 1.5, border = "white", name = NULL) {
  grid::grid.circle(x, y, r = grid::unit(diameter / 2, "mm"), default.units = "native", name = name,
    gp = grid::gpar(fill = fill, col = border, lwd = .35))
}
box_mm <- function(box, w, h) box * c(w, h, w, h)
lookup_colour <- function(value, palette) {
  value <- other(value); value[!value %in% names(palette)] <- "Other"
  unname(palette[value])
}
year_group <- function(x) {
  v <- suppressWarnings(as.integer(x)); out <- rep("Other", length(v))
  ranges <- list(c(2000, 2004), c(2005, 2008), c(2009, 2012), c(2013, 2016), c(2017, 2020), c(2021, 2024))
  labs <- c("2000–2004", "2005–2008", "2009–2012", "2013–2016", "2017–2020", "2021–2024")
  for (i in seq_along(ranges)) out[which(v >= ranges[[i]][1] & v <= ranges[[i]][2])] <- labs[i]
  out
}
amr_group <- function(x) {
  v <- suppressWarnings(as.integer(x)); out <- rep("Other", length(v))
  for (pair in list(list(0, 0, "0"), list(1, 3, "1–3"), list(4, 6, "4–6"), list(7, 10, "7–10"), list(11, 13, "11–13"))) out[which(v >= pair[[1]] & v <= pair[[2]])] <- pair[[3]]
  out
}
metadata_for <- function(d, ids) {
  fields <- c("isolate_id", "run_accession", "sample_accession", "country", "year", "source_group", "lineage", "sequence_type", "AMRFinder_class_count", "XDR_genomic_profile", "plot_label")
  out <- as.data.frame(matrix(NA_character_, length(ids), length(fields)), stringsAsFactors = FALSE)
  names(out) <- fields; out$isolate_id <- ids
  for (i in seq_along(ids)) {
    source <- if (ids[i] %in% rownames(d$meta)) d$meta else d$extras_meta
    if (!is.null(source) && ids[i] %in% rownames(source)) for (f in intersect(fields, names(source))) out[i, f] <- as.character(source[ids[i], f])
  }
  out$isolate_id <- ids; out$plot_label[blank(out$plot_label)] <- ids[blank(out$plot_label)]
  rownames(out) <- ids; out
}
track_values <- function(d, ids, track, cfg) {
  m <- metadata_for(d, ids)
  if (track == "Country") { v <- other(m$country); pal <- cfg$palette$country; v[!v %in% names(pal)] <- "Other" }
  else if (track == "Source") { v <- other(m$source_group); pal <- cfg$palette$source }
  else if (track == "Year") { v <- year_group(m$year); pal <- cfg$palette$year }
  else if (track == "Lineage") { v <- other(m$lineage); pal <- cfg$palette$lineage }
  else if (track == "ST") { v <- other(m$sequence_type); pal <- cfg$palette$st }
  else if (track == "AMR classes") { v <- amr_group(m$AMRFinder_class_count); pal <- cfg$palette$amr_count }
  else if (track %in% c("XDR profile", "CARB-R marker", "MAC-R marker")) {
    col <- c(`XDR profile` = "XDR_genomic_profile", `CARB-R marker` = "CARB_R_genomic_marker", `MAC-R marker` = "MAC_R_genomic_marker")[[track]]
    v <- rep(NA_character_, length(ids)); ii <- match(ids, rownames(d$genomic_flags)); valid <- !is.na(ii)
    v[valid] <- as.character(d$genomic_flags[ii[valid], col])
    v <- ifelse(v == "1", "Profile / marker", ifelse(v == "0", "Other profile", "Other")); v[is.na(v)] <- "Other"
    pal <- cfg$palette$status
  } else if (track %in% ISANGI_AGENTS) {
    v <- rep("Other", length(ids)); ii <- match(ids, rownames(d$phenotype)); ok <- !is.na(ii)
    v[ok] <- other(d$phenotype[ii[ok], track]); pal <- cfg$palette$phenotype
  } else {
    assert(track %in% colnames(d$genes), paste("Unknown track:", track))
    v <- rep("Other", length(ids)); ii <- match(ids, rownames(d$genes)); ok <- !is.na(ii)
    v[ok] <- ifelse(d$genes[ii[ok], track] == 1, "Detected", "Not detected"); pal <- cfg$palette$gene
  }
  list(values = v, palette = pal, colours = lookup_colour(v, pal))
}
legend_block <- function(title, palette, values = NULL, marker = FALSE) {
  if (!is.null(values)) palette <- palette[names(palette) %in% unique(other(values))]
  list(title = title, palette = palette, marker = marker)
}
wrap_mm <- function(x, width, size, cfg) {

  chars <- max(9L, floor(width / (size * cfg$font_scale * .19)))
  paste(strwrap(x, width = chars), collapse = "\n")
}
draw_legends <- function(blocks, x, top, width, bottom, cfg, fontsize = cfg$fonts$legend) {
  blocks <- blocks[vapply(blocks, function(b) length(b$palette) > 0, logical(1))]
  pitch <- (fontsize * cfg$font_scale * .3528) * 1.38
  heading_pitch <- (fontsize * cfg$font_scale * .3528) * 1.60
  calculate <- function() {
    total <- 0
    for (b in blocks) {
      hd <- wrap_mm(b$title, width, fontsize + .3, cfg)
      total <- total + heading_pitch * (1L + lengths(regmatches(hd, gregexpr("\n", hd, fixed = TRUE))))
      for (lab in names(b$palette)) {
        s <- wrap_mm(lab, width - 5, fontsize, cfg)
        total <- total + pitch * (1L + lengths(regmatches(s, gregexpr("\n", s, fixed = TRUE))))
      }; total <- total + pitch * .55
    }; total
  }
  required <- calculate()
  assert(required <= top - bottom + .01, sprintf("Legend needs %.1f mm but has %.1f mm. Increase height/legend width or reduce cfg$fonts$legend.", required, top-bottom))
  y <- top
  for (b in blocks) {
    hd <- wrap_mm(b$title, width, fontsize + .3, cfg)
    text_mm(hd, x, y, cfg, size = fontsize + .3, just = c("left", "top"), bold = TRUE)
    y <- y - heading_pitch * (1 + lengths(regmatches(hd, gregexpr("\n", hd, fixed = TRUE))))
    for (lab in names(b$palette)) {
      s <- wrap_mm(lab, width - 5, fontsize, cfg)
      if (b$marker) point_mm(x + 1.25, y - pitch*.28, b$palette[[lab]], diameter = 1.55, border = "#59606A")
      else rect_mm(x + 1.25, y - pitch*.28, 2.5, 2.4, b$palette[[lab]], "#BDC2C7", .25)
      text_mm(s, x + 4.8, y + .3, cfg, size = fontsize, just = c("left", "top"))
      y <- y - pitch * (1 + lengths(regmatches(s, gregexpr("\n", s, fixed = TRUE))))
    }; y <- y - pitch*.55
  }; invisible(y)
}
metadata_legends <- function(d, ids, tracks, cfg) {
  lapply(tracks, function(k) { t <- track_values(d, ids, k, cfg); legend_block(if (k == "XDR profile") "XDR genomic profile" else k, t$palette, t$values) })
}
row_y <- function(n, b) b[2] + b[4] - (seq_len(n) - .5) * b[4] / n
nice_scale <- function(x) {
  power <- 10^floor(log10(max(x, 1e-12)))
  z <- c(1, 2, 5, 10) * power; z[which.min(abs(z-x))]
}
draw_tree <- function(v, box, d, cfg, tip_diameter = cfg$tree$detail_tip_diameter_mm, tip_track = "Country", scale = TRUE) {
  n <- length(v$ids); yy <- function(y) box[2] + box[4] - (y - .5) * box[4] / n
  xx <- function(x) box[1] + x / (v$max_x * 1.035) * box[3]
  s <- v$segments
  if (nrow(s)) line_mm(xx(s[,1]), yy(s[,2]), xx(s[,3]), yy(s[,4]), cfg$palette$ink, cfg$tree$branch_lwd, "tree-branches")
  colours <- track_values(d, v$ids, tip_track, cfg)$colours
  point_mm(xx(v$tip_x), yy(v$tip_y), colours, tip_diameter, name = "tree-tip-symbols")
  if (scale) {
    sc <- nice_scale(v$max_x * .25); sx <- xx(0)
    line_mm(sx, box[2]-3, xx(sc), box[2]-3, cfg$palette$ink, .8)
    text_mm(format(sc, scientific = FALSE, trim = TRUE), (sx+xx(sc))/2, box[2]-6, cfg, size = 7.0, just = "centre")
  }
}
draw_heatmap <- function(d, ids, columns, box, cfg, label_size = cfg$fonts$heat_column, cell_text = "none") {
  nr <- length(ids); nc <- length(columns); xx <- box[1] + (seq_len(nc)-.5)*box[3]/nc; yy <- row_y(nr, box)
  for (j in seq_along(columns)) {
    co <- track_values(d, ids, columns[j], cfg)$colours
    rect_mm(rep(xx[j], nr), yy, box[3]/nc*.94, box[4]/nr * if (nr <= 45) .94 else 1.0, co, name = paste0("heat-", make.names(columns[j])))
    if (cell_text != "none" && columns[j] %in% ISANGI_AGENTS) {
      lab <- if (cell_text == "zone") as.character(d$zones[ids, columns[j]]) else as.character(d$phenotype[ids, columns[j]])
      ink <- ifelse(d$phenotype[ids, columns[j]] == "R", "white", cfg$palette$ink)
      text_mm(lab, xx[j], yy, cfg, size = cfg$fonts$cell, colour = ink, just = "centre", name = paste0("cell-text-", columns[j]))
    }
  }
  text_mm(columns, xx, box[2]-2.5, cfg, size = label_size, just = "right", rot = 90, name = "heatmap-column-labels")
  invisible(list(x = xx, y = yy))
}
isangi_figure <- function(name, area, draw, width, height, caption = "", mode = "native_R_rebuild", source_data = list(), note = "") {
  structure(list(name = name, area = area, draw = draw, width = width, height = height,
    caption = caption, mode = mode, source_data = source_data, note = note), class = "isangi_figure")
}
plot.isangi_figure <- function(x, ...) {
  grid::grid.newpage(); x$draw(x$width * 25.4, x$height * 25.4); invisible(x)
}
print.isangi_figure <- function(x, ...) plot.isangi_figure(x, ...)
finish_figure <- function(d, cfg, name, draw, width = 10.4, height = 7.1, area = NULL,
                          source_data = list(), mode = "native_R_rebuild", note = "") {
  r <- d$registry[d$registry$name == name, , drop = FALSE]
  if (nrow(r)) { width <- r$width_inches[1]; height <- r$height_inches[1]; area <- area %||% r$area[1] }
  if (!is.null(cfg$sizes[[name]])) { width <- cfg$sizes[[name]][1]; height <- cfg$sizes[[name]][2] }
  caption <- d$captions[[name]] %||% ""
  original_caption <- caption
  if (grepl("^Figure_1Bi", name)) {
    caption <- gsub("CAAP2A, the first blood-culture isolate in the analysis", cfg$timeline$reference, caption, fixed = TRUE)
    caption <- gsub("CAAP2A", cfg$timeline$reference, caption, fixed = TRUE)
    caption <- paste(caption, "Coincident collection-month/source/SNP-distance observations are combined; point area increases with the number of contributing isolates.")
  }
  if (name %in% c("Figure_2_phenotypic_susceptibility", "Figure_S14_comparator_phenotypes") && cfg$phenotype$cell_text != "zone") {
    caption <- gsub("Cells show inhibition-zone diameters (mm)", if (cfg$phenotype$cell_text == "category") "Cell text shows the recorded susceptibility category" else "Cells have no numeric annotation", caption, fixed = TRUE)
    caption <- paste(caption, "Cell annotation mode:", cfg$phenotype$cell_text, "; all recorded zones remain in the accompanying numeric table.")
  }
  if (grepl("^Figure_S(9_|10_|11_|12_|13_)", name)) {
    caption <- gsub("log(1+distance)", switch(cfg$matrix$transform, log1p="log(1+distance)", sqrt="square-root", linear="linear"), caption, fixed=TRUE)
  }
  if (!cfg$panels$use_frozen_membership || cfg$genes$genes_per_block != 28L) {
    if (grepl("^Figure_S(5_|6_|7_|8_|12_|13_)", name)) {
      caption <- paste0("DISPLAY PARTITION CHANGED: use source_data/panel_summary.csv and panel_isolate_index.csv for current membership and revise the associated document caption. Archived caption for reference: ", caption)
      note <- paste(note, "The current display partition differs from the frozen document; the accompanying source index is authoritative.")
    }
  }
  f <- isangi_figure(name, area %||% "supplementary", draw, width, height, caption, mode, source_data, note)
  if (!identical(caption, original_caption)) f$original_caption <- original_caption
  f
}

build_circular <- function(d, cfg) {
  v <- layout_tree(d$cc_tree); st <- v$structure; n <- length(v$ids)
  draw <- function(w, h) {
    native_page(w, h, cfg); on.exit(grid::popViewport())
    gap <- cfg$circular$gap_degrees; assert(gap >= 6 && gap <= 90, "Circular gap must be 6–90 degrees")
    step <- (360-gap)/n*pi/180; sign <- if (cfg$circular$clockwise) -1 else 1
    angles <- rep(NA_real_, st$nn)
    angles[v$nodes] <- pi/2 + sign * (gap*pi/360 + (seq_len(n)-.5)*step)
    aa <- function(node) {
      if (!is.na(angles[node])) return(angles[node])
      ch <- st$children[[as.character(node)]]; z <- vapply(ch, aa, numeric(1))
      angles[node] <<- (z[1]+tail(z,1))/2; angles[node]
    }; aa(st$root)
    r <- cfg$circular$inner_radius + (cfg$circular$tree_radius-cfg$circular$inner_radius) * st$depth / max(st$depth)
    cx <- cfg$circular$centre_fraction[1]*w; cy <- cfg$circular$centre_fraction[2]*h
    mm <- min(h*cfg$circular$radius_fraction, w*.348)/1.12
    px <- function(rad, angle) cx + mm*rad*cos(angle); py <- function(rad, angle) cy + mm*rad*sin(angle)
    sector <- function(r0, r1, a0, a1, fill, border = NA, lwd = .15) {
      a <- seq(a0, a1, length.out = max(4L, ceiling(abs(a1-a0)*120)))
      poly_mm(c(px(r1,a),rev(px(r0,a))), c(py(r1,a),rev(py(r0,a))), fill, border, lwd)
    }
    li <- other(d$meta[v$ids, "lineage"]); rr <- rle(li); starts <- cumsum(c(1L, head(rr$lengths,-1L)))
    for (i in seq_along(starts)) {
      lo <- starts[i]; hi <- lo + rr$lengths[i]-1
      a0 <- angles[v$nodes[lo]]-sign*step/2; a1 <- angles[v$nodes[hi]]+sign*step/2
      sector(.092, cfg$circular$shade_radius, a0, a1, grDevices::adjustcolor(lookup_colour(rr$values[i],cfg$palette$lineage),alpha.f=.38))
    }
    for (node in unique(st$tree$edge[,1])) {
      ch <- st$children[[as.character(node)]]
      a <- seq(angles[ch[1]], angles[tail(ch,1)], length.out = max(4L, ceiling(abs(angles[ch[1]]-angles[tail(ch,1)])*120)))
      grid::grid.lines(px(r[node],a), py(r[node],a), default.units="native", gp=grid::gpar(col=cfg$palette$ink,lwd=cfg$circular$branch_lwd))
      for (u in ch) line_mm(px(r[node],angles[u]),py(r[node],angles[u]),px(r[u],angles[u]),py(r[u],angles[u]),cfg$palette$ink,cfg$circular$branch_lwd)
    }
    point_mm(px(r[v$nodes],angles[v$nodes]), py(r[v$nodes],angles[v$nodes]),
      lookup_colour(d$meta[v$ids,"sequence_type"],cfg$palette$st),cfg$circular$tip_diameter_mm,"#45494E","sequence-type-tips")
    tracks <- c("Country","Source","Year","AMR classes")
    assert(length(cfg$circular$ring_radii)==4 && length(cfg$circular$ring_labels)==4, "Four circular annotation rings are required")
    for (j in seq_along(tracks)) {
      rad <- cfg$circular$ring_radii[j]; co <- track_values(d,v$ids,tracks[j],cfg)$colours
      for (i in seq_len(n)) sector(rad-cfg$circular$ring_width/2,rad+cfg$circular$ring_width/2,
        angles[v$nodes[i]]-step*.49,angles[v$nodes[i]]+step*.49,co[i],"white",.10)
      text_mm(cfg$circular$ring_labels[j],cx,cy+rad*mm,cfg,size=cfg$circular$ring_label_font,just="centre",name=paste0("ring-label-",j))
    }
    blocks <- c(list(legend_block("Sequence type",cfg$palette$st,other(d$meta[v$ids,"sequence_type"]),TRUE),
      legend_block("Lineage",cfg$palette$lineage,li)),metadata_legends(d,v$ids,tracks,cfg))
    draw_legends(blocks,w*.776,h*.973,w*.208,h*.035,cfg,fontsize=7.2)
    sc <- cfg$circular$scale_length; length_mm <- mm*(cfg$circular$tree_radius-cfg$circular$inner_radius)*sc/max(st$depth)
    line_mm(w*.09,h*.035,w*.09+length_mm,h*.035,cfg$palette$ink,.8)
    text_mm(format(sc,trim=TRUE),w*.09+length_mm/2,h*.018,cfg,size=7.2,just="centre")
  }
  finish_figure(d,cfg,"Figure_3A_CC25_circular",draw,source_data=list(tree_tip_order=data.frame(position=seq_len(n),isolate_id=v$ids)))
}

build_st335_main <- function(d,cfg) {
  v <- layout_tree(d$st_tree); columns <- c(cfg$tree$main_tracks,d$main_genes,"XDR profile","CARB-R marker","MAC-R marker")
  draw <- function(w,h) {
    native_page(w,h,cfg); on.exit(grid::popViewport())
    draw_tree(v,box_mm(cfg$tree$main_tree_box,w,h),d,cfg,cfg$tree$main_tip_diameter_mm)
    draw_heatmap(d,v$ids,columns,box_mm(cfg$tree$main_heat_box,w,h),cfg,label_size=cfg$fonts$heat_column)
    blocks <- c(metadata_legends(d,v$ids,c("Country","Year"),cfg),list(
      legend_block("Genomic determinant",cfg$palette$gene,c("Detected","Not detected")),
      legend_block("Genomic status",cfg$palette$status,c("Profile / marker","Other profile"))))
    lb <- box_mm(cfg$tree$main_legend_box,w,h)
    draw_legends(blocks,lb[1],lb[2]+lb[4],lb[3],lb[2],cfg)
  }
  finish_figure(d,cfg,"Figure_4A_ST335_genomic_heatmap",draw,source_data=list(heatmap=figure_heatmap_table(d,v$ids,columns,cfg)))
}
figure_heatmap_table <- function(d,ids,cols,cfg) {
  out <- data.frame(isolate_id=ids,run_accession=metadata_for(d,ids)$run_accession,check.names=FALSE)
  for (col in cols) out[[col]] <- track_values(d,ids,col,cfg)$values
  out
}
build_matched <- function(d,cfg) {
  ids <- d$phenotype$isolate_id[d$phenotype$group=="Isangi"]; v <- layout_tree(d$st_tree,ids)
  draw <- function(w,h) {
    native_page(w,h,cfg);on.exit(grid::popViewport())
    b <- box_mm(c(.025,.23,.195,.70),w,h);draw_tree(v,b,d,cfg,cfg$tree$matched_tip_diameter_mm)
    text_mm(d$meta[v$ids,"plot_label"],w*.232,row_y(length(v$ids),b),cfg,size=cfg$fonts$tip,name="isolate-accession-labels")
    draw_heatmap(d,v$ids,ISANGI_AGENTS,box_mm(c(.525,.23,.342,.70),w,h),cfg)
    draw_legends(list(legend_block("Recorded AST",cfg$palette$phenotype,unlist(d$phenotype[v$ids,ISANGI_AGENTS],use.names=FALSE))),w*.892,h*.89,w*.098,h*.36,cfg)
  }
  finish_figure(d,cfg,"Figure_4B_ST335_matched_phenotypes",draw,source_data=list(heatmap=figure_heatmap_table(d,v$ids,ISANGI_AGENTS,cfg)))
}
build_phenotype <- function(d,cfg,group="Isangi") {
  ids <- d$phenotype$isolate_id[d$phenotype$group==group]
  name <- if(group=="Isangi")"Figure_2_phenotypic_susceptibility" else "Figure_S14_comparator_phenotypes"
  draw <- function(w,h) {
    native_page(w,h,cfg);on.exit(grid::popViewport())
    b <- box_mm(c(.23,.24,.63,.70),w,h)
    draw_heatmap(d,ids,ISANGI_AGENTS,b,cfg,label_size=cfg$fonts$heat_column,cell_text=cfg$phenotype$cell_text)
    text_mm(ids,b[1]-2.5,row_y(length(ids),b),cfg,size=8,just="right",name="phenotype-isolate-labels")
    draw_legends(list(legend_block("Recorded AST",cfg$palette$phenotype,unlist(d$phenotype[ids,ISANGI_AGENTS],use.names=FALSE))),w*.888,h*.90,w*.104,h*.54,cfg,fontsize=8)
    if(cfg$phenotype$cell_text=="zone")text_mm("Cell values:\nzone diameter\n(mm)",w*.888,h*.53,cfg,size=8,just=c("left","top"))
  }
  finish_figure(d,cfg,name,draw,source_data=list(recorded_calls=d$phenotype[ids,,drop=FALSE],recorded_zones=d$zones[ids,,drop=FALSE]))
}
build_panel_index <- function(d,cfg,pop) {
  tr <- if(pop=="CC25")d$cc_tree else d$st_tree;v<-layout_tree(tr);ps<-d$panel_summary[d$panel_summary$population==pop,,drop=FALSE]
  name <- if(pop=="CC25")"Figure_S5_CC25_panel_index" else "Figure_S7_ST335_panel_index"
  draw <- function(w,h) {
    native_page(w,h,cfg);on.exit(grid::popViewport());b<-box_mm(c(.025,.10,.50,.85),w,h)
    draw_tree(v,b,d,cfg,tip_diameter=1.0)
    n<-length(v$ids);yy<-function(p)b[2]+b[4]-(p-.5)*b[4]/n
    for(i in seq_len(nrow(ps))){q<-ps[i,];lo<-as.integer(q$first_tip_position);hi<-as.integer(q$last_tip_position)
      rect_mm(w*.556,(yy(lo)+yy(hi))/2,w*.012,(hi-lo+1)*b[4]/n,lookup_colour(q$lineage,cfg$palette$lineage),"white",.3)
      text_mm(sprintf("%s  |  lineage %s  |  n = %s",q$panel_id,q$lineage,q$n_isolates),w*.577,(yy(lo)+yy(hi))/2,cfg,size=8.4)
      line_mm(w*.548,yy(hi+.5),w*.97,yy(hi+.5),cfg$palette$grid,.35)
    }
  };finish_figure(d,cfg,name,draw,source_data=list(panel_index=ps))
}
build_tree_detail <- function(d,cfg,pop,pid,block=1L) {
  ids<-d$panels$isolate_id[d$panels$panel_id==pid];tr<-if(pop=="CC25")d$cc_tree else d$st_tree;v<-layout_tree(tr,ids)
  if(pop=="CC25"){cols<-cfg$tree$detail_cc25_tracks;name<-paste0("Figure_S6_",pid,"_labelled_tree")}
  else{
    chunks<-split(d$supp_genes,ceiling(seq_along(d$supp_genes)/cfg$genes$genes_per_block))
    assert(block<=length(chunks),"Invalid determinant block");cols<-c("Year","Source",chunks[[block]])
    name<-paste0("Figure_S8_",pid,"_genes_",block)
  }
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport());b<-box_mm(c(.022,.205,.203,.745),w,h)
    draw_tree(v,b,d,cfg);text_mm(d$meta[v$ids,"plot_label"],w*.235,row_y(length(v$ids),b),cfg,size=cfg$fonts$tip,name="isolate-accession-labels")
    if(pop=="CC25"){
      draw_heatmap(d,v$ids,cols,box_mm(c(.523,.205,.20,.745),w,h),cfg,label_size=cfg$fonts$heat_column)
      bl<-metadata_legends(d,v$ids,c("Country","Source","Year","Lineage","AMR classes","XDR profile"),cfg)
      draw_legends(bl,w*.756,h*.95,w*.227,h*.10,cfg)
    }else{
      draw_heatmap(d,v$ids,cols,box_mm(c(.521,.205,.353,.745),w,h),cfg,label_size=cfg$fonts$heat_column)
      bl<-c(metadata_legends(d,v$ids,c("Country","Year","Source"),cfg),list(legend_block("Determinant",cfg$palette$gene,c("Detected","Not detected"))))
      draw_legends(bl,w*.900,h*.95,w*.095,h*.10,cfg,fontsize=6.8)
    }
  }
  finish_figure(d,cfg,name,draw,source_data=list(heatmap=figure_heatmap_table(d,v$ids,cols,cfg)))
}

snp_transform <- function(x, kind) switch(kind,log1p=log1p(x),sqrt=sqrt(x),linear=x,stop("Invalid matrix transform"))
build_snp <- function(d,cfg,pop,pid=NULL,complete=FALSE) {
  detail<-!is.null(pid)
  if(detail){ids<-d$panels$isolate_id[d$panels$panel_id==pid];parent<-if(pop=="CC25")d$cc else d$st;D<-parent[ids,ids,drop=FALSE];name<-paste0(if(pop=="CC25")"Figure_S12_" else "Figure_S13_",pid,"_SNP")}
  else if(complete){assert(!is.null(d$full),"Complete CC25 numeric matrix was not supplied");D<-d$full;name<-"Figure_S9_CC25_complete_331"}
  else if(pop=="CC25"){D<-d$cc;name<-"Figure_S10_CC25_tree_aligned_327"}
  else{D<-d$st;name<-"Figure_S11_ST335_tree_aligned_224"}
  ids<-rownames(D);n<-length(ids);high<-if(pop=="CC25")cfg$matrix$max_cc25 else cfg$matrix$max_st335
  assert(max(D)<=high,"Matrix contains values above the configured colour maximum; increase it explicitly")
  pal<-grDevices::colorRampPalette(cfg$palette$snp)(256)
  scaled<-snp_transform(D,cfg$matrix$transform)/snp_transform(high,cfg$matrix$transform)
  co<-matrix(pal[pmax(1L,pmin(256L,1L+floor(scaled*255)))],n,n)
  tracks<-cfg$matrix$metadata;md<-metadata_for(d,ids)
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport())
    left<-if(detail)w*.28 else w*.095;bottom<-if(detail)h*.205 else h*.16
    side<-min(w*(if(detail).52 else .54),h-bottom-18)
    xx<-left+(seq_len(n)-.5)*side/n;yy<-bottom+side-(seq_len(n)-.5)*side/n
    if(cfg$matrix$draw_cells_as_vectors){

      rect_mm(rep(xx,each=n),rep(yy,times=n),side/n,side/n,as.vector(co),name="SNP-distance-cells")
    }else grid::grid.raster(co,x=left+side/2,y=bottom+side/2,width=side,height=side,default.units="native",interpolate=FALSE,name="SNP-distance-raster")
    rect_mm(left+side/2,bottom+side/2,side,side,NA,cfg$palette$ink,.45)
    strip_h<-1.9
    for(j in seq_along(tracks)){
      t<-track_values(d,ids,tracks[j],cfg)
      strip_y<-bottom+side+2+(length(tracks)-j+.5)*strip_h
      rect_mm(xx,rep(strip_y,n),side/n,strip_h*.88,t$colours)
      text_mm(if(tracks[j]=="XDR profile")"XDR" else tracks[j],left-1.8,strip_y,cfg,size=6.7,just="right")
      rect_mm(rep(left+side+2+(j-.5)*1.0,n),yy,0.88,side/n,t$colours)
    }
    if(detail){
      labels<-if(cfg$matrix$detail_labels=="isolate_accession")md$plot_label else ids
      text_mm(labels,left-2,yy,cfg,size=cfg$fonts$matrix_label,just="right",name="matrix-row-labels")
      text_mm(ids,xx,bottom-2.5,cfg,size=cfg$fonts$matrix_label,rot=90,just="right",name="matrix-column-labels")
    }else{
      tick<-unique(c(1L,if(n>=cfg$matrix$overview_tick_every)seq.int(cfg$matrix$overview_tick_every,n,by=cfg$matrix$overview_tick_every) else integer(),n))
      text_mm(as.character(tick),xx[tick],bottom-3,cfg,size=7.2,just="centre")
      text_mm(as.character(tick),left-2.5,yy[tick],cfg,size=7.2,just="right")
      text_mm("Isolate position in accompanying matrix index",left-11,bottom+side/2,cfg,size=8,rot=90,just="centre")
    }

    cbx<-left;cby<-h*.042;cbw<-side;cbh<-2.8
    rect_mm(cbx+(seq_len(256)-.5)*cbw/256,rep(cby,256),cbw/256+.01,cbh,pal)
    ticks<-c(0,1,2,5,10,20,50,100,if(pop=="CC25")c(500,1000,5000,high) else high)
    ticks<-unique(ticks[ticks<=high]);tx<-cbx+snp_transform(ticks,cfg$matrix$transform)/snp_transform(high,cfg$matrix$transform)*cbw
    for(j in seq_along(ticks)){line_mm(tx[j],cby-cbh/2,tx[j],cby-cbh/2-1,cfg$palette$ink,.4)
      text_mm(format(ticks[j],trim=TRUE,scientific=FALSE),tx[j],cby-cbh/2-2.6,cfg,size=6.4,just=if(j==1)"left" else if(j==length(ticks))"right" else "centre")}
    text_mm(paste0("SNP distance (",switch(cfg$matrix$transform,log1p="log scale",sqrt="square-root scale",linear="linear scale"),")"),cbx+cbw/2,cby+4.6,cfg,size=7.3,just="centre")
    draw_legends(metadata_legends(d,ids,tracks,cfg),w*.862,h*.95,w*.128,h*.12,cfg,fontsize=6.7)
  }
  finish_figure(d,cfg,name,draw,source_data=list(matrix_index=data.frame(matrix_position=seq_len(n),md,check.names=FALSE)))
}

axis_frame <- function(box, cfg, x_ticks, x_labels, y_ticks, y_labels, xlim, ylim,
                       xlab = "", ylab = "", rotate_x = 0) {
  X <- function(x) box[1]+(x-xlim[1])/diff(xlim)*box[3]
  Y <- function(y) box[2]+(y-ylim[1])/diff(ylim)*box[4]
  for (j in seq_along(y_ticks)) {
    line_mm(box[1],Y(y_ticks[j]),box[1]+box[3],Y(y_ticks[j]),cfg$palette$grid,.35)
    text_mm(y_labels[j],box[1]-2,Y(y_ticks[j]),cfg,size=cfg$fonts$axis,just="right")
  }
  line_mm(box[1],box[2],box[1]+box[3],box[2],cfg$palette$ink,.65)
  line_mm(box[1],box[2],box[1],box[2]+box[4],cfg$palette$ink,.65)
  text_mm(x_labels,X(x_ticks),box[2]-3,cfg,size=cfg$fonts$axis,rot=rotate_x,
    just=if(rotate_x==0)c("centre","top") else c("right","centre"))
  if(nzchar(xlab))text_mm(xlab,box[1]+box[3]/2,box[2]-if(rotate_x==0)11 else 21,cfg,size=9.0,just="centre")
  if(nzchar(ylab))text_mm(ylab,box[1]-13,box[2]+box[4]/2,cfg,size=9.0,rot=90,just="centre")
  list(X=X,Y=Y)
}
build_bsi <- function(d,cfg) {
  q<-d$annual;years<-sort(unique(as.integer(q$year)))
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport())
    b<-c(w*.11,h*.20,w*.64,h*.74)
    xy<-axis_frame(b,cfg,years,years,seq(0,100,25),paste0(seq(0,100,25),"%"),
      range(years)+c(-.65,.65),c(0,100),"Collection year","Percentage of isolates")
    for(y in years){
      previous<-0
      for(gr in c("XDR","MDR","Other")){
        p<-as.numeric(q$percent[q$year==y&q$resistance_group==gr]);if(!length(p))p<-0
        assert(length(p)==1L,"Duplicate year/category in BSI summary")
        rect_mm(xy$X(y),xy$Y(previous+p/2),b[3]/(length(years)+.3)*.80,p*b[4]/100,
          cfg$palette$bsi[[gr]],"white",.20);previous<-previous+p
      }
    }
    draw_legends(list(legend_block("Recorded profile",cfg$palette$bsi)),w*.79,h*.76,w*.19,h*.30,cfg)
  }
  finish_figure(d,cfg,"Figure_1A_BSI_resistance",draw,area="main",source_data=list(annual_summary=q))
}
build_timeline <- function(d,cfg,zoom=FALSE) {
  q<-d$timeline
  if(zoom)q<-q[substr(q$collection_month,1,4)=="2020",,drop=FALSE]
  assert(nrow(q)>0,"Empty timeline subset")
  groups<-split(seq_len(nrow(q)),paste(q$collection_month,q$source,q$SNP_distance_to_reference,sep="|"))
  pp<-do.call(rbind,lapply(groups,function(ii){z<-q[ii[1],,drop=FALSE];z$n_coincident<-length(ii);z$isolate_ids<-paste(q$isolate_id[ii],collapse=";");z}))
  if(zoom){limits<-as.Date(c("2020-01-01","2020-12-31"));ticks<-seq(as.Date("2020-01-01"),as.Date("2020-12-01"),by="month");labs<-format(ticks,"%b")}
  else {limits<-as.Date(c(paste0(cfg$timeline$first_year,"-01-01"),paste0(cfg$timeline$last_year,"-12-31")));ticks<-seq(limits[1],limits[2],by="6 months");labs<-format(ticks,"%b %Y")}
  ymax<-max(25,ceiling(max(q$SNP_distance_to_reference)/5)*5)
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport())
    b<-c(w*.09,h*.31,w*.875,h*.64)
    xy<-axis_frame(b,cfg,as.numeric(ticks),labs,seq(0,ymax,5),seq(0,ymax,5),as.numeric(limits),c(-.6,ymax+1),
      "Collection month",paste0("SNP distance from ",cfg$timeline$reference),rotate_x=if(zoom)0 else 45)
    colours<-lookup_colour(pp$source,c(cfg$palette$timeline,Other="#BFC3C9"))
    diameter<-sqrt(cfg$timeline$point_min_mm^2+cfg$timeline$point_count_mm2*(pp$n_coincident-1))
    point_mm(xy$X(as.numeric(pp$plot_date)),xy$Y(pp$SNP_distance_to_reference),colours,diameter,border="white")

    levels<-names(cfg$palette$timeline);x<-w*.10;y<-h*.045
    for(j in seq_along(levels)){point_mm(x,y,cfg$palette$timeline[[levels[j]]],2.1);text_mm(levels[j],x+3.2,y,cfg,size=7.4);x<-x+w*.23}
  }
  name<-if(zoom)"Figure_1Bii_timeline_2020" else "Figure_1Bi_timeline"
  pp$plot_date<-NULL
  finish_figure(d,cfg,name,draw,area="main",source_data=list(plotted_coincident_groups=pp),
    note="The 15th of each month is a drawing coordinate only. Coincident month/source/distance observations are combined by point area; no date jitter is added.")
}

kruskal_snp <- function(D) {
  n<-nrow(D);assert(n>=1L,"Empty network matrix")
  if(n==1L)return(data.frame(from=integer(),to=integer(),snp=numeric()))
  ij<-which(upper.tri(D),arr.ind=TRUE);ord<-order(D[ij],ij[,1],ij[,2]);ij<-ij[ord,,drop=FALSE]
  parent<-seq_len(n);rank<-integer(n)
  find<-function(x){while(parent[x]!=x){parent[x]<<-parent[parent[x]];x<-parent[x]};x}
  out<-matrix(0,n-1,3);k<-0L
  for(i in seq_len(nrow(ij))){a<-ij[i,1];b<-ij[i,2];ra<-find(a);rb<-find(b);if(ra==rb)next
    if(rank[ra]<rank[rb])parent[ra]<-rb else {parent[rb]<-ra;if(rank[ra]==rank[rb])rank[ra]<-rank[ra]+1L}
    k<-k+1L;out[k,]<-c(a,b,D[a,b]);if(k==n-1L)break
  }
  assert(k==n-1L,"Could not form a spanning tree");out<-as.data.frame(out);names(out)<-c("from","to","snp");out
}
zero_components <- function(D) {
  n<-nrow(D);seen<-logical(n);groups<-list()
  for(i in seq_len(n))if(!seen[i]){
    queue<-i;seen[i]<-TRUE;group<-integer()
    while(length(queue)){v<-queue[1];queue<-queue[-1];group<-c(group,v);nb<-which(D[v,]==0 & !seen)
      if(length(nb)){seen[nb]<-TRUE;queue<-c(queue,nb)}}
    groups[[length(groups)+1L]]<-sort(group)
  };groups
}
force_layout_R <- function(n,edges,cfg,radii=rep(.025,n)) {
  if(n==1L)return(matrix(c(.5,.5),1L,dimnames=list(NULL,c("x","y"))))

  old<-if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))get(".Random.seed",envir=.GlobalEnv) else NULL
  on.exit(if(!is.null(old))assign(".Random.seed",old,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(cfg$seed);pos<-matrix(stats::runif(n*2,-.5,.5),n,2);k<-cfg$network$k
  for(iter in seq_len(cfg$network$iterations)){
    dx<-outer(pos[,1],pos[,1],"-");dy<-outer(pos[,2],pos[,2],"-");ds<-sqrt(dx^2+dy^2);diag(ds)<-Inf;ds<-pmax(ds,1e-6)
    f<-cbind(rowSums(dx*k*k/(ds*ds)),rowSums(dy*k*k/(ds*ds)))
    if(nrow(edges))for(e in seq_len(nrow(edges))){a<-edges$from[e];b<-edges$to[e];delta<-pos[a,]-pos[b,];dist<-sqrt(sum(delta^2))+1e-8

      force<-delta*dist/k/max(edges$snp[e],1);f[a,]<-f[a,]-force;f[b,]<-f[b,]+force}
    norm<-sqrt(rowSums(f*f));step<-.075*(1-iter/(cfg$network$iterations+1));pos<-pos+f/pmax(norm,1e-8)*pmin(norm,step)
    pos<-sweep(pos,2,colMeans(pos),"-")
  }
  pos<-pos/max(diff(range(pos[,1])),diff(range(pos[,2])),1e-8)*.72
  for(iter in seq_len(cfg$network$collision_iterations)){
    changed<-FALSE
    for(a in seq_len(n-1L))for(b in seq.int(a+1L,n)){
      delta<-pos[b,]-pos[a,];dist<-sqrt(sum(delta^2));required<-radii[a]+radii[b]+.012
      if(dist<required){if(dist<1e-8){delta<-c(1e-8,0);dist<-1e-8};push<-delta/dist*(required-dist)*.51
        pos[a,]<-pos[a,]-push;pos[b,]<-pos[b,]+push;changed<-TRUE}}
    if(!changed)break
  }

  cbind(x=pos[,1],y=pos[,2])
}
network_data <- function(d,cfg) {
  D<-d$malawi;ids<-rownames(D)
  assert(cfg$network$collapse%in%c("zero_components","none"),"network$collapse must be zero_components or none")
  groups<-if(cfg$network$collapse=="none")as.list(seq_along(ids)) else zero_components(D)
  reps<-vapply(groups,`[`,integer(1),1L);names_nodes<-sprintf("N%03d",seq_along(groups))
  nodes<-do.call(rbind,lapply(seq_along(groups),function(i){ii<-groups[[i]];m<-D[ii,ii,drop=FALSE]
    data.frame(node_id=names_nodes[i],representative_id=ids[reps[i]],n_isolates=length(ii),within_node_max_SNP=max(m),
      all_pairs_zero=all(m==0),grouping_rule=cfg$network$collapse,stringsAsFactors=FALSE)}))
  membership<-do.call(rbind,lapply(seq_along(groups),function(i){ii<-groups[[i]];data.frame(node_id=names_nodes[i],
    isolate_id=ids[ii],run_accession=d$meta[ids[ii],"run_accession"],source=d$timeline$source[match(ids[ii],d$timeline$isolate_id)],stringsAsFactors=FALSE)}))
  edges<-kruskal_snp(D[reps,reps,drop=FALSE]);edges$from_node<-names_nodes[edges$from];edges$to_node<-names_nodes[edges$to]
  edges$from_isolate<-ids[reps[edges$from]];edges$to_isolate<-ids[reps[edges$to]]
  radii<-cfg$network$base_node_radius*sqrt(nodes$n_isolates)
  pos<-force_layout_R(nrow(nodes),edges,cfg,radii)
  if(!is.null(cfg$network$coordinates)){
    xy<-read_tab(path.expand(cfg$network$coordinates));need_cols(xy,c("node_id","x","y"),"Network coordinates");unique_ids(xy$node_id,"Coordinate nodes")
    assert(setequal(xy$node_id,names_nodes),"Edited network coordinates must include every node exactly once")
    xy<-xy[match(names_nodes,xy$node_id),];pos<-cbind(number(xy$x,"Network x"),number(xy$y,"Network y"))
  }
  nodes$x<-pos[,1];nodes$y<-pos[,2];nodes$radius<-radii
  list(nodes=nodes,edges=edges,membership=membership)
}
build_network <- function(d,cfg) {
  q<-network_data(d,cfg);nodes<-q$nodes;edges<-q$edges
  nonclique<-nodes[!nodes$all_pairs_zero,,drop=FALSE]
  note<-paste0("Undirected minimum spanning tree; deterministic Kruskal ties, representative isolates retained in input order. Layout coordinates are arbitrary and not a SNP scale. Grouping: ",cfg$network$collapse,".")
  if(nrow(nonclique))note<-paste0(note," Zero-distance connected components are not necessarily all-pairs-identical: ",
    paste(paste0(nonclique$node_id," (n=",nonclique$n_isolates,", maximum within-node ",nonclique$within_node_max_SNP," SNP)"),collapse="; "),
    ". Use cfg$network$collapse='none' for an uncollapsed view. Node diagnostics and pairwise values are exported.")
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport())
    b<-c(w*.05,h*.07,w*.63,h*.87);pos<-as.matrix(nodes[,c("x","y")]);r<-nodes$radius
    xmin<-min(pos[,1]-r);xmax<-max(pos[,1]+r);ymin<-min(pos[,2]-r);ymax<-max(pos[,2]+r)
    scale<-min(b[3]/(xmax-xmin+.06),b[4]/(ymax-ymin+.06))
    xx<-b[1]+b[3]/2+(pos[,1]-(xmin+xmax)/2)*scale;yy<-b[2]+b[4]/2+(pos[,2]-(ymin+ymax)/2)*scale;rr<-r*scale
    if(nrow(edges))for(e in seq_len(nrow(edges))){a<-edges$from[e];bb<-edges$to[e];line_mm(xx[a],yy[a],xx[bb],yy[bb],"#828A91",.8)
      if(cfg$network$label_edges){x<-mean(xx[c(a,bb)]);y<-mean(yy[c(a,bb)]);lab<-as.character(edges$snp[e]);rect_mm(x,y,2.3+nchar(lab)*1.0,3.4,"white");text_mm(lab,x,y,cfg,size=7.3,just="centre")}}
    for(i in seq_len(nrow(nodes))){
      v<-q$membership$source[q$membership$node_id==nodes$node_id[i]];tt<-table(factor(v,levels=names(cfg$palette$timeline)));cum<-c(0,cumsum(tt)/sum(tt))*2*pi
      for(j in seq_along(tt))if(tt[j]>0){a<-seq(cum[j],cum[j+1],length.out=max(5L,ceiling(100*tt[j]/sum(tt))));
        poly_mm(c(xx[i],xx[i]+rr[i]*cos(a),xx[i]),c(yy[i],yy[i]+rr[i]*sin(a),yy[i]),cfg$palette$timeline[[names(tt)[j]]],"white",.35)}
      a<-seq(0,2*pi,length.out=120);grid::grid.lines(xx[i]+rr[i]*cos(a),yy[i]+rr[i]*sin(a),default.units="native",gp=grid::gpar(col="#6B747D",lwd=.55))
      if(cfg$network$label_node_counts&&nodes$n_isolates[i]>1L)text_mm(nodes$n_isolates[i],xx[i],yy[i],cfg,size=7.5,just="centre",bold=TRUE)
    }
    draw_legends(list(legend_block("Isolate source",cfg$palette$timeline,marker=TRUE)),w*.73,h*.89,w*.245,h*.47,cfg,fontsize=7.5)
    text_mm("Node area: isolate count\nEdge label: SNP distance\nUndirected genomic links",w*.73,h*.40,cfg,size=7.3,just=c("left","top"))

  }
  f<-finish_figure(d,cfg,"Figure_1C_Malawi_minimum_spanning_tree",draw,area="main",source_data=q,note=note)

  f$original_caption<-f$caption
  f$caption<-paste0("Figure 1C. Minimum spanning tree of the 74 Malawian ST335 isolates. ",
    if(cfg$network$collapse=="none")"Each node is one isolate. " else "Nodes group connected components of zero-SNP links; connected membership does not require every pair to be zero SNPs apart. ",
    "Node area denotes isolate count; sectors indicate specimen source. Edges are weighted by the supplied SNP distance between the recorded representative isolates and are undirected, not transmission events. Node membership, within-node distance ranges and representative IDs are provided with the figure.")
  f
}

legacy_path <- function(cfg,key) {
  p<-cfg$legacy[[key]];assert(!is.null(p)&&length(p)==1L&&nzchar(p),paste("Set cfg$legacy$",key,sep=""))
  p<-path.expand(p);assert(file.exists(p),paste("Legacy input not found:",p));normalizePath(p)
}
legacy_table <- function(cfg,key,columns) {
  p<-legacy_path(cfg,key);x<-read_tab(p);need_cols(x,columns,p);assert(nrow(x)>0,paste("Empty input:",p));x
}
build_retained <- function(d,cfg,name) {
  entry<-d$legacy[[which(vapply(d$legacy,function(a)identical(a$name,name),logical(1)))]]
  assert(!is.null(entry),paste("No artwork record for",name))
  im<-png::readPNG(entry$path,native=FALSE);nr<-dim(im)[1];nc<-dim(im)[2]
  if(cfg$legacy$allow_title_masks && length(entry$title_masks))for(msk in entry$title_masks){
    z<-as.integer(unlist(msk));assert(length(z)==4L&&all(z>=0)&&z[3]<=nc&&z[4]<=nr,"Title mask outside artwork")

    im[seq.int(z[2]+1L,z[4]),seq.int(z[1]+1L,z[3]),]<-1
  }
  raster<-as.raster(im)
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport())
    scale<-min(w/nc,h/nr)
    grid::grid.raster(raster,x=w/2,y=h/2,width=nc*scale,height=nr*scale,default.units="native",interpolate=TRUE,name="retained-source-artwork")
  }
  note<-paste0("Retained original artwork: ",entry$source,"; native ",nc," x ",nr," pixels. No underlying measurements inferred; vector exports contain this raster. ",
    if(length(entry$title_masks)&&cfg$legacy$allow_title_masks)"Only the documented plot-title rectangles were masked as in the previous release." else "No scientific content redrawn.")
  f<-finish_figure(d,cfg,name,draw,width=entry$width_inches,height=entry$width_inches*nr/nc,area=entry$area,
    mode="retained_source_artwork",note=note,source_data=list(artwork_provenance=data.frame(source_path=entry$path,
      native_width_px=nc,native_height_px=nr,md5=unname(tools::md5sum(entry$path)),stringsAsFactors=FALSE)))
  f
}
build_map_native <- function(d,cfg) {
  assert(requireNamespace("sf",quietly=TRUE)&&requireNamespace("ggplot2",quietly=TRUE),"Native map requires optional packages sf and ggplot2")
  fp<-legacy_path(cfg,"map_gpkg");layers<-sf::st_layers(fp)$name
  need<-c("boundary","waterways","sampling_sites","hospital");assert(all(need%in%layers),"Map GPKG requires boundary, waterways, sampling_sites and hospital layers")
  obj<-setNames(lapply(need,function(n){x<-sf::st_read(fp,layer=n,quiet=TRUE);assert(nrow(x)>0,paste("Empty map layer",n))
    assert(!is.na(sf::st_crs(x)),paste("Undeclared CRS for",n));assert(all(sf::st_is_valid(x)),paste("Invalid map geometry in",n))
    sf::st_transform(x,cfg$legacy$map_epsg)}),need)
  need_cols(obj$sampling_sites,c("site_id","category"),"Sampling sites");unique_ids(obj$sampling_sites$site_id,"Map sites")
  need_cols(obj$hospital,"label","Hospital points")
  obj$sampling_sites$category<-other(obj$sampling_sites$category)
  categories<-unique(obj$sampling_sites$category)
  cols<-setNames(grDevices::hcl.colors(length(categories),"Dark 3"),categories)
  if("Salmonella Isangi"%in%categories)cols["Salmonella Isangi"]<-cfg$palette$country[["Malawi"]]
  if("Other"%in%categories)cols["Other"]<-cfg$palette$country[["Other"]]
  links<-NULL
  if(!is.null(cfg$legacy$map_membership)){
    links<-legacy_table(cfg,"map_membership",c("site_id","isolate_id","run_accession"));unique_ids(links$isolate_id,"Mapped isolates")
    assert(all(links$site_id%in%obj$sampling_sites$site_id),"Isolate-site mapping refers to absent sites")
    ii<-match(links$isolate_id,d$meta$isolate_id);ok<-!is.na(ii)
    assert(all(d$meta$run_accession[ii[ok]]==links$run_accession[ok]),"Map accession linkage conflicts with master metadata")
  }
  gp<-ggplot2::ggplot()+ggplot2::geom_sf(data=obj$boundary,fill="#F3F5F2",colour="#AAB3B9",linewidth=.35)+
    ggplot2::geom_sf(data=obj$waterways,colour="#6CA4BD",linewidth=.28)+
    ggplot2::geom_sf(data=obj$sampling_sites,ggplot2::aes(colour=category),size=2.0,alpha=.9)+
    ggplot2::geom_sf(data=obj$hospital,shape=23,fill="#E5AD43",colour="#6B5940",size=3)+
    ggplot2::geom_sf_text(data=obj$hospital,ggplot2::aes(label=label),nudge_y=150,size=2.5)+
    ggplot2::scale_colour_manual(values=cols,name="Sample category")+
    ggplot2::coord_sf(crs=sf::st_crs(cfg$legacy$map_epsg),datum=NA)+
    ggplot2::labs(title=NULL,subtitle=NULL,x=NULL,y=NULL)+
    ggplot2::theme_void(base_size=10,base_family=cfg$font_family)+
    ggplot2::theme(legend.position="right",legend.title=ggplot2::element_text(size=9),legend.text=ggplot2::element_text(size=8),
      plot.margin=ggplot2::margin(8,8,8,8),plot.background=ggplot2::element_rect(fill="white",colour=NA))
  grob<-ggplot2::ggplotGrob(gp);draw<-function(w,h){grid::grid.draw(grob)}
  sd<-list(sampling_sites=sf::st_drop_geometry(obj$sampling_sites),hospital=sf::st_drop_geometry(obj$hospital))
  if(!is.null(links))sd$isolate_site_membership<-links
  f<-finish_figure(d,cfg,"Figure_1D_Blantyre_map",draw,area="main",source_data=sd,
    note="Native R map from the explicitly configured GPKG. No guessed coordinates, CRS or river names; uniform waterway width unless modified by the author. New layout, not pixel-equivalence to archived map.")
  f$editable_plot<-gp;f
}
build_burden_native <- function(d,cfg) {
  q<-legacy_table(cfg,"ncbi_counts",c("isolate_id","group","n_amr_genes","snapshot_date"));unique_ids(q$isolate_id,"NCBI comparison")
  q$n_amr_genes<-number(q$n_amr_genes,"NCBI gene counts");assert(all(q$n_amr_genes>=0&q$n_amr_genes==round(q$n_amr_genes)),"AMR gene counts must be non-negative integers")
  assert(length(unique(q$snapshot_date))==1L&&!any(blank(q$snapshot_date)),"One explicit frozen snapshot date is required")
  q$group<-other(q$group);lev<-unique(q$group);assert(length(lev)>=2L,"AMR comparison requires at least two observed groups")
  colours<-setNames(grDevices::hcl.colors(length(lev),"Set 2"),lev)
  summary<-do.call(rbind,lapply(lev,function(g){v<-q$n_amr_genes[q$group==g];data.frame(group=g,n=length(v),median=stats::median(v),q1=unname(stats::quantile(v,.25)),q3=unname(stats::quantile(v,.75)),snapshot_date=q$snapshot_date[1])}))
  densities<-lapply(lev,function(g){v<-q$n_amr_genes[q$group==g];if(length(unique(v))<2L)return(NULL);stats::density(v,from=0,to=max(v),n=512)})
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport());b<-c(w*.14,h*.23,w*.80,h*.68)
    high<-max(5,ceiling(max(q$n_amr_genes)/5)*5);xy<-axis_frame(b,cfg,seq_along(lev),vapply(lev,wrap_mm,character(1),width=b[3]/length(lev)*.90,size=8,cfg=cfg),
      pretty(c(0,high)),pretty(c(0,high)),c(.4,length(lev)+.6),c(0,high*1.05),"","Number of AMR genes")
    for(i in seq_along(lev)){v<-q$n_amr_genes[q$group==lev[i]];den<-densities[[i]]
      if(!is.null(den)){half<-den$y/max(den$y)*.36;poly_mm(c(xy$X(i-half),rev(xy$X(i+half))),c(xy$Y(den$x),rev(xy$Y(den$x))),colours[i],"#6C747D",.5)}
      else line_mm(xy$X(i-.25),xy$Y(v[1]),xy$X(i+.25),xy$Y(v[1]),colours[i],3)
      line_mm(xy$X(i),xy$Y(summary$q1[i]),xy$X(i),xy$Y(summary$q3[i]),"#3B4148",1.7)
      point_mm(xy$X(i),xy$Y(summary$median[i]),"white",1.8,"#3B4148")
      text_mm(paste0("n = ",format(length(v),big.mark=",")),xy$X(i),b[2]+b[4]+6,cfg,size=7.5,just="centre")
    }
  }
  finish_figure(d,cfg,"Figure_3B_AMR_gene_burden",draw,area="main",source_data=list(comparison_summary=summary),
    note=paste("Native kernel-density comparison from the frozen supplied export dated",q$snapshot_date[1],". Groups and denominators are observed, not copied or inferred from the old image."))
}

validate_plasmids <- function(cfg) {
  p<-legacy_table(cfg,"plasmids",c("plasmid_id","isolate_id","length_bp","replicon","assembly_accession"));unique_ids(p$plasmid_id,"Plasmids")
  p$length_bp<-number(p$length_bp,"Plasmid lengths");assert(all(p$length_bp>0&p$length_bp==round(p$length_bp)),"Invalid plasmid lengths")
  f<-legacy_table(cfg,"plasmid_features",c("plasmid_id","start","end","strand","gene","is_amr","region_id"))
  r<-legacy_table(cfg,"plasmid_regions",c("plasmid_id","region_id","start","end"))
  for(qn in c("f","r")){
    q<-get(qn);q$start<-number(q$start,paste(qn,"start"));q$end<-number(q$end,paste(qn,"end"))
    lim<-p$length_bp[match(q$plasmid_id,p$plasmid_id)]
    assert(all(is.finite(lim))&&all(q$start>=1&q$end>=q$start&q$end<=lim),paste("Invalid 1-based inclusive coordinates in",qn));assign(qn,q)
  }
  assert(all(f$strand%in%c("+","-"))&&all(f$is_amr%in%c("0","1")),"Invalid plasmid gene strand/is_amr")
  unique_ids(paste(r$plasmid_id,r$region_id,sep="|"),"Plasmid regions")
  matched<-!blank(f$region_id);ii<-match(paste(f$plasmid_id[matched],f$region_id[matched],sep="|"),paste(r$plasmid_id,r$region_id,sep="|"))
  assert(!anyNA(ii)&&all(f$start[matched]>=r$start[ii]&f$end[matched]<=r$end[ii]),"Gene coordinates exceed/omit named resistance region")
  a<-NULL
  if(!is.null(cfg$legacy$plasmid_alignments)){
    a<-legacy_table(cfg,"plasmid_alignments",c("query_plasmid","query_start","query_end","target_plasmid","target_start","target_end","percent_identity"))
    for(n in c("query_start","query_end","target_start","target_end","percent_identity"))a[[n]]<-number(a[[n]],n)
    for(side in c("query","target")){lim<-p$length_bp[match(a[[paste0(side,"_plasmid")]],p$plasmid_id)]
      assert(all(is.finite(lim))&&all(a[[paste0(side,"_start")]]>=1&a[[paste0(side,"_end")]]>=1&a[[paste0(side,"_start")]]<=lim&a[[paste0(side,"_end")]]<=lim),"Plasmid alignment endpoint outside molecule")}
    assert(all(a$percent_identity>=0&a$percent_identity<=100),"Alignment identity outside 0–100")
  };list(plasmids=p,features=f,regions=r,alignments=a)
}
gene_arrow <- function(x0,x1,y,strand,fill,label,cfg,half_height=1.8) {
  tip<-min(abs(x1-x0)*.45,2.0)
  if(strand=="+")poly_mm(c(x0,x1-tip,x1,x1-tip,x0),c(y-half_height,y-half_height,y,y+half_height,y+half_height),fill,"#555E67",.25)
  else poly_mm(c(x1,x0+tip,x0,x0+tip,x1),c(y-half_height,y-half_height,y,y+half_height,y+half_height),fill,"#555E67",.25)
  if(!blank(label))text_mm(label,(x0+x1)/2,y+half_height+2.5,cfg,size=6.4,rot=38,just="left")
}
build_plasmids_native <- function(d,cfg) {
  q<-validate_plasmids(cfg);p<-q$plasmids;f<-q$features;r<-q$regions;regions<-unique(r$region_id)
  genes<-unique(other(f$gene));pal<-setNames(grDevices::hcl.colors(length(genes),"Dynamic"),genes)
  draw<-function(w,h){
    native_page(w,h,cfg);on.exit(grid::popViewport())
    left<-w*.18;right<-w*.97;span<-right-left;maxlen<-max(p$length_bp)
    yy<-seq(h*.89,h*.64,length.out=nrow(p));X<-function(bp)left+span*bp/maxlen
    if(!is.null(q$alignments))for(i in seq_len(nrow(q$alignments))){a<-q$alignments[i,];u<-match(a$query_plasmid,p$plasmid_id);v<-match(a$target_plasmid,p$plasmid_id)
      poly_mm(c(X(a$query_start),X(a$query_end),X(a$target_end),X(a$target_start)),c(yy[u],yy[u],yy[v],yy[v]),
        grDevices::adjustcolor("#88ABA0",alpha.f=.15+.35*a$percent_identity/100),NA)}
    for(i in seq_len(nrow(p))){line_mm(X(1),yy[i],X(p$length_bp[i]),yy[i],"#9AABB6",3)
      text_mm(paste0(p$plasmid_id[i],"\n",p$replicon[i]),left-3,yy[i],cfg,size=7.0,just="right")
      ff<-f[f$plasmid_id==p$plasmid_id[i]&f$is_amr=="1",,drop=FALSE]
      for(j in seq_len(nrow(ff)))gene_arrow(X(ff$start[j]),X(ff$end[j]),yy[i],ff$strand[j],pal[[other(ff$gene[j])]],"",cfg,1.4)
    }
    scale_bp<-nice_scale(maxlen/5);line_mm(left,h*.595,X(scale_bp),h*.595,"#30373F",.7)
    text_mm(paste0(format(scale_bp/1000,trim=TRUE)," kb"),mean(c(left,X(scale_bp))),h*.577,cfg,size=7,just="centre")
    ncol<-min(2,length(regions));nrow<-ceiling(length(regions)/ncol);bh<-h*.47/nrow
    for(ri in seq_along(regions)){
      col<-(ri-1)%%ncol;row<-(ri-1)%/%ncol;x0<-w*(.03+col*.49);top<-h*.50-row*bh;bw<-w*.46
      rr<-r[r$region_id==regions[ri],,drop=FALSE];active<-p$plasmid_id[p$plasmid_id%in%rr$plasmid_id]
      text_mm(regions[ri],x0+1,top+3,cfg,size=7.8,bold=TRUE)
      max_region<-max(rr$end-rr$start+1);y0<-seq(top-10,top-bh+8,length.out=length(active))
      for(j in seq_along(active)){
        rrj<-rr[rr$plasmid_id==active[j],];ff<-f[f$plasmid_id==active[j]&!is.na(f$region_id)&f$region_id==regions[ri],,drop=FALSE]
        gx<-function(z)x0+bw*.23+(z-rrj$start+1)/max_region*bw*.75
        text_mm(active[j],x0,y0[j],cfg,size=6.1)
        line_mm(gx(rrj$start),y0[j],gx(rrj$end),y0[j],"#B5BDC3",.6)
        for(k in seq_len(nrow(ff)))gene_arrow(gx(ff$start[k]),gx(ff$end[k]),y0[j],ff$strand[k],pal[[other(ff$gene[k])]],ff$gene[k],cfg,1.1)
      }
    }
  }
  sd<-q[!vapply(q,is.null,logical(1))]
  finish_figure(d,cfg,"Figure_5_plasmid_comparison",draw,area="main",source_data=sd,
    note=paste("Native coordinate-based plasmid reconstruction.",if(is.null(q$alignments))"No alignment input: NO homology ribbons were drawn." else "Ribbons use supplied alignment blocks only.","Visual layout differs from archived comparison; genes and orientation remain source-defined."))
}
build_plasmid_network_native <- function(d,cfg) {
  nodes<-legacy_table(cfg,"plasmid_nodes",c("plasmid_id","replicon"));unique_ids(nodes$plasmid_id,"Plasmid network nodes")
  edges<-legacy_table(cfg,"plasmid_edges",c("from","to","dcj_indel_distance"))
  edges$dcj_indel_distance<-number(edges$dcj_indel_distance,"DCJ-Indel distances");assert(all(edges$dcj_indel_distance>=0),"Negative DCJ-Indel distance")
  ee<-data.frame(from=match(edges$from,nodes$plasmid_id),to=match(edges$to,nodes$plasmid_id),snp=edges$dcj_indel_distance)
  assert(!anyNA(ee)&&all(ee$from!=ee$to),"Invalid network edge endpoints")
  unique_ids(paste(pmin(ee$from,ee$to),pmax(ee$from,ee$to),sep="|"),"Undirected plasmid edges")
  pos<-force_layout_R(nrow(nodes),ee,cfg,rep(.025,nrow(nodes)))
  if(!is.null(cfg$legacy$plasmid_positions)){
    xy<-legacy_table(cfg,"plasmid_positions",c("plasmid_id","x","y"));unique_ids(xy$plasmid_id,"Plasmid positions")
    assert(setequal(xy$plasmid_id,nodes$plasmid_id),"Plasmid positions do not cover every node")
    xy<-xy[match(nodes$plasmid_id,xy$plasmid_id),];pos<-cbind(number(xy$x,"Plasmid x"),number(xy$y,"Plasmid y"))
  }
  nodes$x<-pos[,1];nodes$y<-pos[,2];replicons<-unique(other(nodes$replicon));pal<-setNames(grDevices::hcl.colors(length(replicons),"Dark 3"),replicons)
  draw<-function(w,h){native_page(w,h,cfg);on.exit(grid::popViewport());b<-c(w*.04,h*.04,w*.70,h*.89)
    scale<-min(b[3]/max(diff(range(pos[,1])),.1),b[4]/max(diff(range(pos[,2])),.1))*.85
    xx<-b[1]+b[3]/2+(pos[,1]-mean(range(pos[,1])))*scale;yy<-b[2]+b[4]/2+(pos[,2]-mean(range(pos[,2])))*scale
    line_mm(xx[ee$from],yy[ee$from],xx[ee$to],yy[ee$to],grDevices::adjustcolor("#657781",alpha.f=.3),.35)
    point_mm(xx,yy,lookup_colour(nodes$replicon,c(pal,Other="#BFC3C9")),2.3)
    highlight<-grepl("^CP028197(\\.1)?$",nodes$plasmid_id)
    if(any(highlight))text_mm(nodes$plasmid_id[highlight],xx[highlight]+2,yy[highlight]+3,cfg,size=7.6,bold=TRUE)
    draw_legends(list(legend_block("Replicon",pal,marker=TRUE)),w*.78,h*.89,w*.20,h*.25,cfg,fontsize=7.1)
  }
  finish_figure(d,cfg,"Figure_S1_plasmid_network",draw,source_data=list(nodes=nodes,edges=edges),
    note="Native R plasmid similarity graph from explicitly provided nodes and DCJ-Indel edges. No new edge threshold; geometry is a force layout, not distance-scaled. Full identifiers and edge values are exported.")
}
build_confocal_native <- function(d,cfg) {
  fp<-legacy_path(cfg,"confocal_manifest");q<-legacy_table(cfg,"confocal_manifest",c("strain","image_path"))
  unique_ids(q$strain,"Confocal image panels")
  paths<-vapply(q$image_path,function(p){if(!grepl("^(/|~|[A-Za-z]:)",p))p<-file.path(dirname(fp),p);p<-path.expand(p);assert(file.exists(p),paste("Confocal image absent",p));p},character(1))
  assert(all(grepl("\\.png$",paths,ignore.case=TRUE)),"Confocal adapter accepts lossless PNG panels; convert original images without scientific alteration first")
  ims<-lapply(paths,png::readPNG);q$resolved_path<-paths;q$md5<-unname(tools::md5sum(paths))
  draw<-function(w,h){native_page(w,h,cfg);on.exit(grid::popViewport());nc<-min(2,nrow(q));nr<-ceiling(nrow(q)/nc)
    for(i in seq_len(nrow(q))){col<-(i-1)%%nc;row<-(i-1)%/%nc;bw<-w/nc*.94;bh<-h/nr*.87;x<-(col+.5)*w/nc;y<-h-(row+.5)*h/nr
      im<-ims[[i]];scale<-min(bw/dim(im)[2],bh/dim(im)[1]);grid::grid.raster(im,x=x,y=y,width=dim(im)[2]*scale,height=dim(im)[1]*scale,default.units="native",interpolate=FALSE)
      text_mm(q$strain[i],x,y+bh/2+3,cfg,size=8.5,just="centre")}
  }
  finish_figure(d,cfg,"Figure_S2_confocal_biofilms",draw,source_data=list(panel_manifest=q),
    note="Native R assembly of real confocal image panels; images remain raster. No intensity adjustment, segmentation, invented scale bar or reanalysis of Z-stacks.")
}
build_biofilm_native <- function(d,cfg) {
  q<-legacy_table(cfg,"biofilm_values",c("strain","replicate","metric","component","value"))
  q$value<-number(q$value,"Biofilm values");assert(all(q$value>=0),"Negative biofilm measurement")
  unique_ids(paste(q$strain,q$replicate,q$metric,q$component,sep="|"),"Biological-replicate biofilm values")

  splitq<-split(seq_len(nrow(q)),paste(q$strain,q$metric,q$component,sep="|"))
  ss<-do.call(rbind,lapply(splitq,function(ii){z<-q[ii[1],c("strain","metric","component"),drop=FALSE];z$n<-length(ii);z$mean<-mean(q$value[ii]);z$sd<-if(length(ii)>1L)stats::sd(q$value[ii]) else NA_real_;z}))
  assert(all(ss$n>=2),"Mean +/- SD panels require at least two independent replicates in every group")
  metric<-unique(q$metric);components<-unique(q$component);strains<-unique(q$strain)
  pal<-c(Curli="#C94E55",Cells="#6FAF71",Cellulose="#5384B4");more<-setdiff(components,names(pal));if(length(more))pal<-c(pal,setNames(grDevices::hcl.colors(length(more),"Dark 3"),more))
  ann<-NULL;if(!is.null(cfg$legacy$biofilm_annotations)){
    ann<-legacy_table(cfg,"biofilm_annotations",c("metric","component","strain_a","strain_b","y","label"));ann$y<-number(ann$y,"Annotation y")
    assert(all(ann$metric%in%metric&ann$component%in%components&ann$strain_a%in%strains&ann$strain_b%in%strains),"Invalid explicit biofilm annotation")}
  draw<-function(w,h){native_page(w,h,cfg);on.exit(grid::popViewport());nc<-min(2,length(metric));nr<-ceiling(length(metric)/nc)
    for(mi in seq_along(metric)){col<-(mi-1)%%nc;row<-(mi-1)%/%nc;bw<-w/nc;bh<-h/nr;origin<-c(col*bw,h-(row+1)*bh)
      b<-c(origin[1]+bw*.14,origin[2]+bh*.27,bw*.70,bh*.61);xpos<-expand.grid(strain=strains,component=components,stringsAsFactors=FALSE)
      sub<-ss[ss$metric==metric[mi],];high<-max(sub$mean+sub$sd)*1.12
      if(!is.null(ann)&&any(ann$metric==metric[mi]))high<-max(high,ann$y[ann$metric==metric[mi]]*1.1)
      ticks<-pretty(c(0,high),n=4);ticks<-ticks[ticks>=0&ticks<=high]
      xy<-axis_frame(b,cfg,seq_len(nrow(xpos)),xpos$strain,ticks,ticks,c(.4,nrow(xpos)+.6),c(0,high),"",metric[mi],rotate_x=50)
      for(j in seq_len(nrow(xpos))){a<-sub[sub$strain==xpos$strain[j]&sub$component==xpos$component[j],];if(!nrow(a))next
        rect_mm(xy$X(j),xy$Y(a$mean/2),b[3]/nrow(xpos)*.68,xy$Y(a$mean)-xy$Y(0),pal[[xpos$component[j]]])
        line_mm(xy$X(j),xy$Y(max(0,a$mean-a$sd)),xy$X(j),xy$Y(a$mean+a$sd),"#36414A",.55)
        line_mm(xy$X(j)-1,xy$Y(a$mean+a$sd),xy$X(j)+1,xy$Y(a$mean+a$sd),"#36414A",.55)}
      text_mm(LETTERS[mi],origin[1]+2,origin[2]+bh*.94,cfg,size=11,bold=TRUE)
      if(!is.null(ann))for(k in which(ann$metric==metric[mi])){
        a<-ann[k,];i<-which(xpos$strain==a$strain_a&xpos$component==a$component);j<-which(xpos$strain==a$strain_b&xpos$component==a$component)
        line_mm(xy$X(i),xy$Y(a$y),xy$X(j),xy$Y(a$y),"#36414A",.55);text_mm(a$label,mean(xy$X(c(i,j))),xy$Y(a$y)+2,cfg,size=7,just="centre")}
    }

    x<-w*.12;for(cmp in components){rect_mm(x,h-3,2.3,2.3,pal[[cmp]]);text_mm(cmp,x+3.8,h-3,cfg,size=7.0);x<-x+w*.24}
  }
  sd<-list(independent_replicate_values=q,summary=ss);if(!is.null(ann))sd$explicit_annotations<-ann
  finish_figure(d,cfg,"Figure_S3_biofilm_quantification",draw,source_data=sd,
    note="Native mean +/- SD from supplied independent replicate-level measurements. Technical Z-stacks are not treated as additional biological replicates. No p-values or stars inferred from the archived image; annotations only from an explicitly supplied table.")
}
kaplan_meier_R <- function(time,event) {
  assert(length(time)==length(event)&&all(is.finite(time))&&all(time>=0)&&all(event%in%c(0,1)),"Invalid survival data")
  tt<-sort(unique(time));out<-data.frame(time_days=0,survival=1,n_at_risk=length(time),events=0,censored=0);s<-1
  for(t in tt){nr<-sum(time>=t);ne<-sum(time==t&event==1);nc<-sum(time==t&event==0);s<-s*(1-ne/nr)
    out<-rbind(out,data.frame(time_days=t,survival=s,n_at_risk=nr,events=ne,censored=nc))};out
}
build_survival_native <- function(d,cfg) {
  q<-legacy_table(cfg,"mouse_survival",c("mouse_id","strain","time_days","event"));unique_ids(q$mouse_id,"Mouse IDs")
  q$time_days<-number(q$time_days,"Follow-up time");q$event<-number(q$event,"Survival event")
  lev<-unique(q$strain);curves<-do.call(rbind,lapply(lev,function(g){z<-q[q$strain==g,];s<-kaplan_meier_R(z$time_days,z$event);s$strain<-g;s}))
  pal<-setNames(grDevices::hcl.colors(length(lev),"Dark 3"),lev);if("14028"%in%lev)pal["14028"]<-"#222222"
  draw<-function(w,h){native_page(w,h,cfg);on.exit(grid::popViewport());b<-c(w*.12,h*.20,w*.61,h*.71);end<-max(q$time_days)
    ticks<-pretty(c(0,end),n=5);ticks<-ticks[ticks>=0&ticks<=end]
    xy<-axis_frame(b,cfg,ticks,ticks,seq(0,1,.25),paste0(seq(0,100,25),"%"),c(0,end),c(0,1.025),"Days elapsed","Survival")
    for(i in seq_along(lev)){z<-curves[curves$strain==lev[i],];for(j in seq.int(2L,nrow(z))){line_mm(xy$X(z$time_days[j-1]),xy$Y(z$survival[j-1]),xy$X(z$time_days[j]),xy$Y(z$survival[j-1]),pal[i],1.1)
      line_mm(xy$X(z$time_days[j]),xy$Y(z$survival[j-1]),xy$X(z$time_days[j]),xy$Y(z$survival[j]),pal[i],1.1)}
      cens<-which(z$censored>0);if(length(cens))line_mm(xy$X(z$time_days[cens]),xy$Y(z$survival[cens])-1,xy$X(z$time_days[cens]),xy$Y(z$survival[cens])+1,pal[i],.8)
    };draw_legends(list(legend_block("Strain",pal,marker=TRUE)),w*.77,h*.87,w*.21,h*.28,cfg)
  }
  finish_figure(d,cfg,"Figure_S4_mouse_survival",draw,source_data=list(mouse_records=q,kaplan_meier_estimates=curves),
    note="Native Kaplan-Meier steps from actual per-animal follow-up and event indicators, with censoring marks. No event times inferred from the published image and no hypothesis tests invented.")
}
legacy_specification <- function() list(
  Figure_1D_Blantyre_map=list(keys="map_gpkg",builder=build_map_native),
  Figure_3B_AMR_gene_burden=list(keys="ncbi_counts",builder=build_burden_native),
  Figure_5_plasmid_comparison=list(keys=c("plasmids","plasmid_features","plasmid_regions"),builder=build_plasmids_native),
  Figure_S1_plasmid_network=list(keys=c("plasmid_nodes","plasmid_edges"),builder=build_plasmid_network_native),
  Figure_S2_confocal_biofilms=list(keys="confocal_manifest",builder=build_confocal_native),
  Figure_S3_biofilm_quantification=list(keys="biofilm_values",builder=build_biofilm_native),
  Figure_S4_mouse_survival=list(keys="mouse_survival",builder=build_survival_native)
)
build_legacy <- function(d,cfg,name) {
  spec<-legacy_specification()[[name]];mode<-cfg$legacy$mode
  assert(mode%in%c("retain","prefer_data","require_data"),"legacy$mode must be retain, prefer_data or require_data")
  available<-all(vapply(spec$keys,function(k)!is.null(cfg$legacy[[k]])&&length(cfg$legacy[[k]])==1&&nzchar(cfg$legacy[[k]]),logical(1)))
  if(mode=="retain" || (mode=="prefer_data"&&!available))return(build_retained(d,cfg,name))
  assert(available,paste("Native reconstruction requires cfg$legacy fields:",paste(spec$keys,collapse=", ")))

  spec$builder(d,cfg)
}

figure_plan <- function(d,cfg) {
  plan<-list()
  add<-function(name,builder,tags){plan[[name]]<<-list(builder=builder,tags=unique(c(tags,if(grepl("^Figure_S",name))"supplementary" else "main")))}
  add("Figure_1A_BSI_resistance",build_bsi,"epidemiology")
  add("Figure_1Bi_timeline",function(d,cfg)build_timeline(d,cfg,FALSE),"epidemiology")
  add("Figure_1Bii_timeline_2020",function(d,cfg)build_timeline(d,cfg,TRUE),"epidemiology")
  add("Figure_1C_Malawi_minimum_spanning_tree",build_network,"epidemiology")
  add("Figure_2_phenotypic_susceptibility",function(d,cfg)build_phenotype(d,cfg,"Isangi"),"phenotype")
  add("Figure_3A_CC25_circular",build_circular,c("trees","cc25"))
  add("Figure_4A_ST335_genomic_heatmap",build_st335_main,c("trees","st335"))
  add("Figure_4B_ST335_matched_phenotypes",build_matched,c("trees","st335","phenotype"))
  for(pop in c("CC25","ST335"))local({p<-pop
    name<-if(p=="CC25")"Figure_S5_CC25_panel_index" else "Figure_S7_ST335_panel_index"
    add(name,function(d,cfg)build_panel_index(d,cfg,p),c("trees",tolower(p)))
  })
  for(pid in unique(d$panels$panel_id))local({id<-pid;pop<-d$panels$population[match(id,d$panels$panel_id)]
    if(pop=="CC25")add(paste0("Figure_S6_",id,"_labelled_tree"),function(d,cfg)build_tree_detail(d,cfg,pop,id),c("trees","cc25"))
    else for(block in seq_len(ceiling(length(d$supp_genes)/cfg$genes$genes_per_block)))local({bb<-block
      add(paste0("Figure_S8_",id,"_genes_",bb),function(d,cfg)build_tree_detail(d,cfg,pop,id,bb),c("trees","st335"))})
    name<-paste0(if(pop=="CC25")"Figure_S12_" else "Figure_S13_",id,"_SNP")
    add(name,function(d,cfg)build_snp(d,cfg,pop,id),c("snp",tolower(pop)))
  })
  if(!is.null(d$full))add("Figure_S9_CC25_complete_331",function(d,cfg)build_snp(d,cfg,"CC25",complete=TRUE),c("snp","cc25"))
  add("Figure_S10_CC25_tree_aligned_327",function(d,cfg)build_snp(d,cfg,"CC25"),c("snp","cc25"))
  add("Figure_S11_ST335_tree_aligned_224",function(d,cfg)build_snp(d,cfg,"ST335"),c("snp","st335"))
  add("Figure_S14_comparator_phenotypes",function(d,cfg)build_phenotype(d,cfg,"Comparator"),"phenotype")
  for(name in names(legacy_specification()))local({id<-name;add(id,function(d,cfg)build_legacy(d,cfg,id),"legacy")})

  ord<-c(d$registry$name[d$registry$name%in%names(plan)],setdiff(names(plan),d$registry$name));plan<-plan[ord]
  for(name in intersect(names(cfg$custom_builders),names(plan)))plan[[name]]$builder<-cfg$custom_builders[[name]]
  if(cfg$check_frozen_results&&cfg$panels$use_frozen_membership&&cfg$genes$genes_per_block==28L){
    expected<-d$registry$name;if(is.null(d$full))expected<-setdiff(expected,"Figure_S9_CC25_complete_331")
    assert(setequal(names(plan),expected),"The R build plan does not cover the frozen embedded figure registry")
  }
  plan
}

build_isangi_figure <- function(name,d,cfg=isangi_config()) {
  p<-figure_plan(d,cfg);assert(name%in%names(p),paste("No builder for",name))
  f<-p[[name]]$builder(d,cfg);assert(inherits(f,"isangi_figure"),"A custom builder must return an isangi_figure")
  f
}
open_figure_device <- function(path,format,width,height,cfg) {
  switch(format,
    tiff=ragg::agg_tiff(path,width=width,height=height,units="in",res=cfg$dpi,compression="lzw",background="white"),
    png=ragg::agg_png(path,width=width,height=height,units="in",res=cfg$preview_dpi,background="white"),
    pdf={if(capabilities("cairo"))grDevices::cairo_pdf(path,width=width,height=height,family=cfg$font_family,onefile=FALSE,bg="white")
      else grDevices::pdf(path,width=width,height=height,family="Helvetica",onefile=FALSE,useDingbats=FALSE,bg="white")},
    svg=svglite::svglite(path,width=width,height=height,bg="white"),
    stop("Unsupported format: ",format,call.=FALSE))
}

save_isangi_figure <- function(f,cfg,out) {
  assert(inherits(f,"isangi_figure"),"Expected an isangi_figure")
  assert(all(is.finite(c(f$width,f$height)))&&f$width>0&&f$height>0,"Invalid physical figure dimensions")
  result<-list();i<-0L;written<-character();complete<-FALSE
  on.exit(if(!complete&&length(written))unlink(written),add=TRUE)
  for(fmt in cfg$formats){
    area<-if(fmt%in%c("pdf","svg"))"vector" else if(fmt=="png")"previews" else f$area
    dir.create(file.path(out,area),recursive=TRUE,showWarnings=FALSE)
    path<-file.path(out,area,paste0(f$name,".",fmt));assert(!file.exists(path),paste("Output already exists:",path))
    tryCatch({
      open_figure_device(path,fmt,f$width,f$height,cfg);dev<-grDevices::dev.cur()
      tryCatch({grid::grid.newpage();f$draw(f$width*25.4,f$height*25.4)},finally={if(dev%in%grDevices::dev.list())grDevices::dev.off(dev)})
    },error=function(e){if(file.exists(path))unlink(path);stop(e)})
    written<-c(written,path)
    i<-i+1L;result[[i]]<-data.frame(figure_id=f$name,mode=f$mode,format=fmt,relative_path=file.path(area,basename(path)),
      width_inches=f$width,height_inches=f$height,dpi=if(fmt=="tiff")cfg$dpi else if(fmt=="png")cfg$preview_dpi else NA,
      bytes=file.info(path)$size,md5=unname(tools::md5sum(path)),stringsAsFactors=FALSE)
  }
  if(cfg$save_grobs){
    g<-grid::grid.grabExpr({grid::grid.newpage();f$draw(f$width*25.4,f$height*25.4)},wrap=TRUE,width=f$width,height=f$height)
    dir.create(file.path(out,"objects"),recursive=TRUE,showWarnings=FALSE)
    objpath<-file.path(out,"objects",paste0(f$name,"_grob.rds"));assert(!file.exists(objpath),paste("Object already exists:",objpath))
    saveRDS(g,objpath,compress="xz");written<-c(written,objpath)
    if(!is.null(f$editable_plot)){ggpath<-file.path(out,"objects",paste0(f$name,"_ggplot.rds"));assert(!file.exists(ggpath),paste("Object already exists:",ggpath));saveRDS(f$editable_plot,ggpath,compress="xz");written<-c(written,ggpath)}
  }
  dir.create(file.path(out,"captions"),recursive=TRUE,showWarnings=FALSE)
  cap_path<-file.path(out,"captions",paste0(f$name,".txt"));assert(!file.exists(cap_path),paste("Caption already exists:",cap_path))
  writeLines(c(f$caption,"",paste0("Technical reproduction note: ",f$note)),cap_path,useBytes=TRUE);written<-c(written,cap_path)
  if(!is.null(f$original_caption)){cp<-file.path(out,"captions",paste0(f$name,"_archived_caption.txt"));assert(!file.exists(cp),paste("Caption already exists:",cp));writeLines(f$original_caption,cp,useBytes=TRUE);written<-c(written,cp)}
  for(n in names(f$source_data))if(is.data.frame(f$source_data[[n]])||is.matrix(f$source_data[[n]])){
    dir.create(file.path(out,"source_data"),recursive=TRUE,showWarnings=FALSE)
    dp<-file.path(out,"source_data",paste0(f$name,"__",n,".csv"));assert(!file.exists(dp),paste("Source table already exists:",dp))
    write_tab(as.data.frame(f$source_data[[n]],check.names=FALSE),dp);written<-c(written,dp)
  }
  complete<-TRUE;do.call(rbind,result)
}
assemble_figures <- function(figures,cfg,names,output_name) {
  assert(all(names%in%base::names(figures)),paste("Assembly requires",paste(names,collapse=", ")))
  parts<-figures[names];width<-max(vapply(parts,`[[`,numeric(1),"width"));gap<-.15
  heights<-vapply(parts,function(f)f$height*width/f$width,numeric(1));height<-sum(heights)+gap*(length(parts)-1)
  draw<-function(w,h){
    top<-1
    for(i in seq_along(parts)){
      hh<-heights[i]/height
      grid::pushViewport(grid::viewport(x=.5,y=top-hh/2,width=1,height=hh))

      parts[[i]]$draw(w,hh*h);grid::popViewport();top<-top-hh-gap/height
    }
  }
  if(!is.null(cfg$sizes[[output_name]])){width<-cfg$sizes[[output_name]][1];height<-cfg$sizes[[output_name]][2]}
  isangi_figure(output_name,"assemblies",draw,width,height,paste(vapply(parts,`[[`,character(1),"caption"),collapse="\n\n"),
    mode=if(any(vapply(parts,function(p)p$mode=="retained_source_artwork",logical(1))))"native_R_assembly_with_retained_artwork" else "native_R_assembly",
    note="Production assembly. Individual plates remain separately editable and are the ones embedded in the clean documents.")
}
escape_html <- function(x) {x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE);x<-gsub(">","&gt;",x,fixed=TRUE);gsub('"',"&quot;",x,fixed=TRUE)}
write_gallery <- function(status,out) {
  cards<-vapply(seq_len(nrow(status)),function(i){s<-status[i,];p<-paste0("previews/",s$figure_id,".png")
    paste0("<section><h2>",escape_html(s$figure_id),"</h2><p class='mode'>",escape_html(s$mode)," | ",escape_html(s$status),"</p>",
      if(file.exists(file.path(out,p)))paste0("<a href='",p,"'><img loading='lazy' src='",p,"'></a>") else "",
      "<p>",escape_html(s$note),"</p></section>")},character(1))
  writeLines(c("<!doctype html><meta charset='utf-8'><title>Isangi native R figure gallery</title>",
    "<style>body{font:16px system-ui;margin:32px;max-width:1250px}section{border-top:1px solid #ccc;margin:36px 0}img{max-width:100%;height:auto}.mode{color:#555}h1{font-size:26px}h2{font-size:19px}</style>",
    "<h1>Isangi figures — native R render</h1><p>Click a preview for full resolution. Vector PDF/SVG files are in vector/. Captions and source tables are separate. Retained artwork is explicitly identified.</p>",cards),file.path(out,"index.html"),useBytes=TRUE)
}
validate_config <- function(cfg) {
  assert(length(cfg$targets)>0L&&all(!blank(cfg$targets)),"Choose at least one figure target")
  assert(all(cfg$formats%in%c("tiff","pdf","svg","png"))&&!anyDuplicated(cfg$formats)&&length(cfg$formats)>0L,"Invalid formats")
  assert(is.numeric(cfg$dpi)&&length(cfg$dpi)==1L&&is.finite(cfg$dpi)&&cfg$dpi>=72,"Invalid TIFF dpi")
  assert(cfg$panels$max_rows>=2&&cfg$genes$genes_per_block>=1,"Invalid panel partition settings")
  assert(cfg$phenotype$cell_text%in%c("zone","category","none"),"cell_text must be zone, category or none")
  assert(cfg$matrix$transform%in%c("log1p","sqrt","linear"),"Invalid SNP colour transform")
  assert(cfg$font_scale>0&&cfg$circular$gap_degrees>=6,"Invalid font/gap controls")
  invisible(cfg)
}
run_isangi <- function(cfg=isangi_config(),check_only=FALSE) {
  validate_config(cfg);check_dependencies(if(check_only)character() else cfg$formats)
  project<-normalizePath(path.expand(cfg$project_root),mustWork=FALSE)
  if(!dir.exists(project))dir.create(project,recursive=TRUE,showWarnings=FALSE)
  assert(dir.exists(project),paste("Could not create project directory:",project))
  cfg$project_root<-project
  if(is.null(cfg$output_dir)){
    parent<-if(grepl("^(/|~|[A-Za-z]:)",cfg$output_parent))path.expand(cfg$output_parent) else file.path(project,cfg$output_parent)
    out<-file.path(parent,paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_nativeR_",Sys.getpid()))
  }else out<-path.expand(cfg$output_dir)
  assert(!dir.exists(out),paste("Refusing to overwrite an existing run:",out))
  dir.create(file.path(out,"audit"),recursive=TRUE,showWarnings=FALSE)

  utils::capture.output(dput(cfg),file=file.path(out,"audit","configuration.R"))
  utils::capture.output(utils::sessionInfo(),file=file.path(out,"audit","sessionInfo_start.txt"))
  src<-file.path(ISANGI_HOME,"build_isangi_figures.R");if(file.exists(src))file.copy(src,file.path(out,"audit","build_isangi_figures.R"))
  message("Input root: ",cfg$input_root,"\nNew output folder: ",out)
  d<-load_isangi(cfg,out);plan<-figure_plan(d,cfg)
  plan_table<-data.frame(figure_id=names(plan),groups=vapply(plan,function(p)paste(p$tags,collapse=";"),character(1)),stringsAsFactors=FALSE)
  write_tab(plan_table,file.path(out,"audit","complete_figure_plan.csv"))
  network_check<-network_data(d,cfg);write_tab(network_check$nodes,file.path(out,"audit","Malawi_zero_component_diagnostics.csv"))
  if(any(!network_check$nodes$all_pairs_zero))message("NOTE: one or more zero-distance connected components are not all-pairs-identical. See audit/Malawi_zero_component_diagnostics.csv and the revised Figure 1C caption.")
  if(check_only){message("PASS: input validation and figure-plan coverage. No figures rendered.");return(invisible(list(output_dir=out,data=d,plan=plan_table,figures=list())))}
  valid_targets<-unique(c("all",names(plan),unlist(lapply(plan,`[[`,"tags"),use.names=FALSE)))
  assert(all(cfg$targets%in%valid_targets),paste("Unrecognised target(s):",paste(setdiff(cfg$targets,valid_targets),collapse=", ")))
  selected<-names(plan)[vapply(names(plan),function(n)"all"%in%cfg$targets||n%in%cfg$targets||any(plan[[n]]$tags%in%cfg$targets),logical(1))]
  figures<-list();status<-list();manifest<-list();idx<-0L
  write_progress<-function(){
    if(length(status))write_tab(do.call(rbind,status),file.path(out,"audit","figure_status.csv"))
    if(length(manifest))write_tab(do.call(rbind,manifest),file.path(out,"audit","output_manifest.csv"))
  }
  for(n in selected){
    message("Building ",n)
    warnings_seen<-character()
    res<-tryCatch(withCallingHandlers({
      f<-build_isangi_figure(n,d,cfg);mf<-save_isangi_figure(f,cfg,out);list(figure=f,manifest=mf)
    },warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
    idx<-idx+1L
    if(inherits(res,"error")){
      status[[idx]]<-data.frame(figure_id=n,status="FAILED",mode="not_generated",note=conditionMessage(res),warnings=paste(warnings_seen,collapse=" | "),stringsAsFactors=FALSE)
      message("FAILED: ",conditionMessage(res))
    }else{
      figures[[n]]<-res$figure;manifest[[length(manifest)+1L]]<-res$manifest
      status[[idx]]<-data.frame(figure_id=n,status="PASS",mode=res$figure$mode,note=res$figure$note,warnings=paste(warnings_seen,collapse=" | "),stringsAsFactors=FALSE)
    };write_progress()
  }
  if(cfg$create_assemblies){
    assemblies<-list(Figure_1_assembled=c("Figure_1A_BSI_resistance","Figure_1Bi_timeline","Figure_1Bii_timeline_2020","Figure_1C_Malawi_minimum_spanning_tree","Figure_1D_Blantyre_map"),
      Figure_3_assembled=c("Figure_3A_CC25_circular","Figure_3B_AMR_gene_burden"),Figure_4_assembled=c("Figure_4A_ST335_genomic_heatmap","Figure_4B_ST335_matched_phenotypes"))
    for(n in names(assemblies))if(all(assemblies[[n]]%in%names(figures))){
      message("Assembling ",n);res<-tryCatch({f<-assemble_figures(figures,cfg,assemblies[[n]],n);mf<-save_isangi_figure(f,cfg,out);list(figure=f,manifest=mf)},error=function(e)e)
      idx<-idx+1L
      if(inherits(res,"error"))status[[idx]]<-data.frame(figure_id=n,status="FAILED",mode="not_generated",note=conditionMessage(res),warnings="",stringsAsFactors=FALSE)
      else {figures[[n]]<-res$figure;manifest[[length(manifest)+1L]]<-res$manifest;status[[idx]]<-data.frame(figure_id=n,status="PASS",mode=res$figure$mode,note=res$figure$note,warnings="",stringsAsFactors=FALSE)}
      write_progress()
    }
  }
  st<-do.call(rbind,status);write_gallery(st,out)
  utils::capture.output(utils::sessionInfo(),file=file.path(out,"audit","sessionInfo_end.txt"))

  lp<-unlist(cfg$legacy[setdiff(names(cfg$legacy),c("mode","map_epsg","allow_title_masks"))],use.names=TRUE)
  lp<-lp[vapply(lp,function(p)file.exists(path.expand(p)),logical(1))]
  if(length(lp))write_tab(data.frame(role=names(lp),path=path.expand(lp),md5=unname(tools::md5sum(path.expand(lp)))),file.path(out,"audit","optional_source_inputs.csv"))
  failed<-sum(st$status=="FAILED");message("\n",sum(st$status=="PASS")," successful outputs; ",failed," failed.\n",file.path(out,"index.html"))
  result<-list(output_dir=out,data=d,figures=if(cfg$return_figures)figures else NULL,status=st,manifest=if(length(manifest))do.call(rbind,manifest) else NULL,configuration=cfg)
  if(failed&&cfg$fail_at_end)stop(failed," figure(s) failed. Successful plates and the exact errors are saved in ",out,". No failed figure is marked complete.",call.=FALSE)
  invisible(result)
}

isangi_self_test <- function(input_root=NULL) {
  tests<-list();check<-function(name,condition){assert(condition,paste("SELF-TEST FAILED:",name));tests[[name]]<<-TRUE}
  cfg<-isangi_config()
  check("limited_identifier_normalisation",identical(clean_id(c("CIV13XE8_S20_L001","44341_3#96","ABC_assembly_contigs.fasta")),c("CIV13XE8","44341_3_96","ABC")))
  check("browser_suffix_normalisation",identical(canonical_filename(c("a(2).CSV","a (1)(1).csv")),c("a.csv","a.csv")))
  check("Other_label",identical(other(c(NA,"Unknown","","Malawi")),c("Other","Other","Other","Malawi")))
  check("year_2000_retained",identical(year_group(c(2000,2004,2005,2024,NA)),c("2000–2004","2000–2004","2005–2008","2021–2024","Other")))
  D<-matrix(c(0,2,2,0),2,dimnames=list(c("a","b"),c("a","b")))
  check("matrix_column_reordering",identical(validate_matrix(D[,2:1,drop=FALSE]),D))
  bad<-D;bad[1,2]<-3;check("asymmetric_matrix_rejected",inherits(try(validate_matrix(bad),silent=TRUE),"try-error"))
  bad<-D;bad[1,2]<-NA;check("missing_distance_not_zero",inherits(try(validate_matrix(bad),silent=TRUE),"try-error"))
  bad<-D;diag(bad)<-1;check("nonzero_diagonal_rejected",inherits(try(validate_matrix(bad),silent=TRUE),"try-error"))
  check("normalisation_collision_rejected",inherits(try(unique_ids(clean_id(c("a#b","a_b")),"test"),silent=TRUE),"try-error"))
  p<-as.data.frame(matrix("S",4,length(ISANGI_AGENTS)),stringsAsFactors=FALSE);names(p)<-ISANGI_AGENTS;p$isolate_id<-paste0("test",1:4);p$group<-"Isangi"
  p[1,c("AMP10","C30","SXT25","PEF5","CPD10")]<-"R"
  p[2,c("AMP10","C30","SXT25","PEF5","CTX5")]<-"R"
  p[3,c("AMP10","C30","SXT25","CIP5","CPD10")]<-"R"
  p[4,]<-p[1,];p$isolate_id[4]<-"test4";p$PEF5[4]<-NA
  ph<-classify_phenotypes(p,cfg)
  check("recorded_phenotype_rule",identical(ph$XDR_phenotype,c(1L,1L,0L,NA_integer_)))
  check("exact_17_agents_no_AMX",length(ISANGI_AGENTS)==17L&&!"AMX"%in%ISANGI_AGENTS&&"CTX5"%in%ISANGI_AGENTS)
  dd<-matrix(c(0,0,1,0,0,0,1,0,0),3,dimnames=list(letters[1:3],letters[1:3]))
  cc<-zero_components(dd);check("zero_connected_not_identical",length(cc)==1L&&max(dd[cc[[1]],cc[[1]]])==1)
  ee<-kruskal_snp(dd);check("true_zero_edges_preserved",nrow(ee)==2L&&sum(ee$snp)==0)
  km<-kaplan_meier_R(c(2,4,4,5),c(1,0,1,0))
  check("Kaplan_Meier_censoring",isTRUE(all.equal(km$survival,c(1,.75,.50,.50))))
  if(requireNamespace("ape",quietly=TRUE)){
    t<-tree_structure(ape::read.tree(text="((a:1,b:2):3,c:4);"));v<-layout_tree(t)
    check("frozen_tree_order",identical(t$ids,c("a","b","c")))
    check("tree_tip_distances",isTRUE(all.equal(as.numeric(v$tip_x),c(4,5,4))))
    vv<-layout_tree(t,c("a","b"));check("subset_preserves_branch_lengths",isTRUE(all.equal(as.numeric(vv$tip_x),c(1,2))))
  }
  if(!is.null(input_root)){
    check_dependencies(character());cfg$input_root<-input_root;cfg$project_root<-tempdir();d<-load_isangi(cfg)
    pl<-figure_plan(d,cfg);check("all_embedded_plates_have_builders",length(pl)==66L)
    check("all_curated_genes_in_supplement",length(d$supp_genes)==55L)
    check("12_CC25_7_ST335_panels",length(unique(d$panels$panel_id[d$panels$population=="CC25"]))==12L&&length(unique(d$panels$panel_id[d$panels$population=="ST335"]))==7L)

    for(n in names(pl)){f<-build_isangi_figure(n,d,cfg);check(paste0("construct_",n),inherits(f,"isangi_figure"))}
  }
  result<-data.frame(test=names(tests),result="PASS",stringsAsFactors=FALSE);print(result,row.names=FALSE);message(length(tests)," self-tests passed.");invisible(result)
}

isangi_cli <- function(args=commandArgs(trailingOnly=TRUE)) {
  help<-c("Native R Isangi figure build", "",
    "Rscript build_isangi_figures.R --project '/path/to/S. Isangi'",
    "Rscript build_isangi_figures.R --project '/path/to/S. Isangi' --check-only",
    "Rscript build_isangi_figures.R --self-test",
    "",
    "Options: --project PATH; --inputs PATH; --source-zip ZIP (repeatable);",
    "         --target all|main|supplementary|trees|cc25|st335|snp|phenotype|epidemiology|legacy|FIGURE_ID",
    "         --formats tiff,pdf,svg,png; --dpi 600; --check-only; --self-test;",
    "         --no-assemblies; --no-grobs; --help",
    "",
    "For editable configuration, source this file in RStudio and call run_isangi(cfg).",
    "No Python, ggtree or Quarto dependency. Rendering never installs packages.")
  if("--help"%in%args||"-h"%in%args){cat(paste(help,collapse="\n"),"\n");return(invisible(NULL))}
  values<-list();flags<-character();i<-1L
  valued<-c("project","inputs","source-zip","target","formats","dpi")
  while(i<=length(args)){
    a<-args[i]
    if(grepl("^--[^=]+=",a)){key<-sub("^--([^=]+)=.*$","\\1",a);value<-sub("^--[^=]+=","",a)}
    else if(sub("^--","",a)%in%valued){key<-sub("^--","",a);i<-i+1L;assert(i<=length(args),paste("Missing value for",a));value<-args[i]}
    else {assert(a%in%c("--check-only","--self-test","--no-assemblies","--no-grobs"),paste("Unrecognised option",a));flags<-c(flags,a);i<-i+1L;next}
    assert(key%in%valued,paste("Unrecognised option",key));values[[key]]<-c(values[[key]],value);i<-i+1L
  }
  cfg<-isangi_config(project_root=tail(values$project%||%getwd(),1L),input_root=tail(values$inputs%||%file.path(ISANGI_HOME,"inputs"),1L))
  if(!is.null(values[["source-zip"]]))cfg$source_zip<-values[["source-zip"]]
  if(!is.null(values$target))cfg$targets<-unlist(strsplit(values$target,",",fixed=TRUE))
  if(!is.null(values$formats))cfg$formats<-strsplit(tail(values$formats,1L),",",fixed=TRUE)[[1]]
  if(!is.null(values$dpi))cfg$dpi<-as.integer(tail(values$dpi,1L))
  cfg$create_assemblies<-!"--no-assemblies"%in%flags;cfg$save_grobs<-!"--no-grobs"%in%flags
  if("--self-test"%in%flags)return(isangi_self_test(cfg$input_root))
  run_isangi(cfg,check_only="--check-only"%in%flags)
}
