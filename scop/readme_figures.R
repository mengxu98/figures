devtools::document()

fig_dir <- "../figures/scop"
fig_dir <- normalizePath(fig_dir)

examples_dpi <- 210

# svglite reports 72 px per inch, so a raster_quality factor of dpi / 72
# renders heatmap rasters at that dpi (600 px for a 2 in panel at 300 dpi).
raster_dpi <- 300
ht_raster_params <- list(raster_quality = raster_dpi / 72)

# Rasterising now lives in thisplot (>= 0.4.6.9000) as rasterise_plot().

for (pkg in c("svglite", "thisplot", "png")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    thisutils::log_message(
      "Package {.pkg {pkg}} is required to write the SVG and dark-mode figures.",
      message_type = "error"
    )
  }
}

save_fig <- function(
  p,
  filename,
  width = 9,
  height = 3,
  units = "in",
  dpi = 210
) {
  svg <- file.path(fig_dir, filename)
  ggsave(
    svg,
    plot = p,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )
  invert_args <- list(path = svg, quiet = TRUE)
  if ("invert_raster" %in% names(formals(thisplot::invert_svg))) {
    invert_args[["invert_raster"]] <- TRUE
  } else if (grepl("<image", readChar(svg, file.info(svg)$size))) {
    thisutils::log_message(
      "  The loaded {.pkg thisplot} cannot invert embedded rasters; {.file {basename(svg)}} keeps its light colours in the dark variant.",
      message_type = "warning"
    )
  }
  dark <- do.call(thisplot::invert_svg, invert_args)
  thisutils::log_message(
    "  Saved: {.file {basename(svg)}} and {.file {basename(dark)}}",
    message_type = "success"
  )
}

external_figures <- c("SCExplorer-1")

save_dark_png <- function(figure) {
  light <- file.path(fig_dir, paste0(figure, ".png"))
  if (!file.exists(light)) {
    thisutils::log_message(
      "  {.file {basename(light)}} is missing; produce that figure outside this script first.",
      message_type = "warning"
    )
    return(invisible(NULL))
  }

  dark <- thisplot::invert_png(light, quiet = TRUE)
  thisutils::log_message(
    "  Saved: {.file {basename(light)}} and {.file {basename(dark)}}",
    message_type = "success"
  )
  invisible(dark)
}

just_margin <- FALSE


thisutils::log_message("=== Section: Data exploration ===")

data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(
  pancreas_sub
)
print(pancreas_sub)

p <- CellDimPlot(
  pancreas_sub,
  group.by = c("CellType", "SubCellType"),
  reduction = "UMAP",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "EDA-1.svg", width = 9, height = 3)

p <- CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  stat.by = "Phase",
  reduction = "UMAP",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "EDA-2.svg", width = 9, height = 3)

p <- FeatureDimPlot(
  pancreas_sub,
  features = c("Sox9", "Neurog3", "Fev", "Rbp4"),
  reduction = "UMAP",
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "EDA-3.svg", width = 9, height = 4)

p <- FeatureDimPlot(
  pancreas_sub,
  features = c("Ins1", "Gcg", "Sst", "Ghrl"),
  compare_features = TRUE,
  label = TRUE,
  label_insitu = TRUE,
  reduction = "UMAP",
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "EDA-4.svg", width = 9, height = 3)

ht <- GroupHeatmap(
  pancreas_sub,
  features = c(
    "Sox9",
    "Anxa2", # Ductal
    "Neurog3",
    "Hes6", # EPs
    "Fev",
    "Neurod1", # Pre-endocrine
    "Rbp4",
    "Pyy", # Endocrine
    "Ins1",
    "Gcg",
    "Sst",
    "Ghrl" # Beta, Alpha, Delta, Epsilon
  ),
  exp_legend_title = "Z-score",
  group.by = c("CellType", "SubCellType"),
  heatmap_palette = "Spectral",
  group_palette = "ChineseSet8",
  cell_annotation = c("Phase", "G2M_score", "Cdh2"),
  cell_annotation_palette = c("Dark2", "ChineseSet8", "ChineseSet8"),
  show_row_names = TRUE,
  nlabel = 0,
  row_names_side = "left",
  add_dot = TRUE,
  dot_size = grid::unit(6, "mm"),
  add_reticle = TRUE,
  width = 2.2,
  height = 2.2
)
if (just_margin) {
  p <- ht$plot &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
} else {
  p <- ht$plot
}
save_fig(
  p,
  "EDA-5.svg",
  width = 8.65,
  height = 4.37
)


thisutils::log_message("=== Section: CellQC ===")

pancreas_sub <- RunCellQC(pancreas_sub)

p <- CellDimPlot(
  pancreas_sub,
  group.by = "CellQC",
  palette = "ChineseSet8",
  reduction = "UMAP"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunCellQC-1.svg", width = 9, height = 3)

p <- CellStatPlot(
  pancreas_sub,
  stat.by = "CellQC",
  group.by = "CellType",
  palette = "ChineseSet8"
) +
  ggplot2::theme(aspect.ratio = 1 / 2)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunCellQC-2.svg", width = 9, height = 3)

p <- CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "db_qc",
    "outlier_qc",
    "umi_qc",
    "gene_qc",
    "mito_qc",
    "ribo_qc",
    "ribo_mito_ratio_qc",
    "species_qc"
  ),
  plot_type = "upset",
  stat_level = "Fail"
) +
  ggplot2::theme(aspect.ratio = 1 / 2)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunCellQC-3.svg", width = 9, height = 3)


thisutils::log_message("=== Section: Integration pipeline ===")

data(panc8_sub)
panc8_sub <- integration_scop(
  srt_merge = panc8_sub,
  batch = "tech",
  integration_method = "Harmony"
)

p <- CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "HarmonyUMAP",
  theme_use = "theme_blank",
  palette = "ChineseSet8"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "integration_scop-1.svg", width = 9, height = 3)


thisutils::log_message("=== Section: Cell projection ===")

genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub[["RNA"]]),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_rename <- RenameFeatures(
  panc8_sub,
  newnames = genenames,
  assays = "RNA"
)
srt_query <- RunKNNMap(
  srt_query = pancreas_sub,
  srt_ref = panc8_rename,
  ref_umap = "HarmonyUMAP2D"
)

p <- ProjectionPlot(
  srt_query = srt_query,
  srt_ref = panc8_rename,
  query_param = list(
    palette = "ChineseSet8",
    cells.highlight = TRUE
  ),
  ref_param = list(
    palette = "ChineseSet8",
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  ),
  query_group = "SubCellType",
  ref_group = "celltype"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunKNNMap-1.svg", width = 7.5, height = 3.2)


thisutils::log_message("=== Section: Cell annotation (bulk) ===")

data("ref_scMCA")
pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  bulk_ref = ref_scMCA,
  filter_lowfreq = 20
)

p <- CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  reduction = "UMAP",
  label = TRUE,
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunKNNPredict-bulk-1.svg", width = 7.5, height = 3)


thisutils::log_message("=== Section: Cell annotation (single-cell) ===")

pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  srt_ref = panc8_rename,
  ref_group = "celltype",
  filter_lowfreq = 20
)

p <- CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  reduction = "UMAP",
  label = TRUE,
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunKNNPredict-scrna-1.svg", width = 7.5, height = 3)

ht <- CellCorHeatmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_rename,
  query_group = "SubCellType",
  ref_group = "celltype",
  nlabel = 3,
  label_by = "row",
  query_group_palette = "ChineseSet8",
  show_row_names = TRUE,
  show_column_names = TRUE,
  width = 4,
  height = 2.5
)
if (just_margin) {
  p <- ht$plot &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
} else {
  p <- ht$plot
}
save_fig(p, "RunKNNPredict-scrna-3.svg", width = 8.7, height = 4.6)


thisutils::log_message("=== Section: CytoTRACE 2 ===")

pancreas_sub <- RunCytoTRACE(
  pancreas_sub,
  species = "Mus_musculus"
)

p <- CytoTRACEPlot(
  pancreas_sub,
  group.by = "CellType",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
p <- p &
  ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )
save_fig(
  p,
  "RunCytoTRACE.svg",
  width = 12.5,
  height = 6
)


thisutils::log_message("=== Section: SCVELO ===")

pancreas_sub <- RunSCVELO(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP"
)

p <- VelocityPlot(
  pancreas_sub,
  reduction = "UMAP",
  group.by = "SubCellType",
  group_palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunSCVELO-1.svg", width = 9, height = 3)

p <- VelocityPlot(
  pancreas_sub,
  reduction = "UMAP",
  plot_type = "stream",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunSCVELO-2.svg", width = 9, height = 3)


thisutils::log_message("=== Section: PAGA analysis ===")

pancreas_sub <- RunPAGA(
  pancreas_sub,
  group.by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP"
)

p <- PAGAPlot(
  pancreas_sub,
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE,
  node_palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunPAGA-1.svg", width = 9, height = 3)


thisutils::log_message("=== Section: Slingshot ===")

pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  start = "Ductal",
  end = c("Alpha", "Beta"),
  show_plot = FALSE
)

p <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  lineages = paste0("Lineage", 1:2),
  reduction = "UMAP",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunSlingshot-1.svg", width = 9, height = 3)

p <- FeatureDimPlot(
  pancreas_sub,
  features = paste0("Lineage", 1:2),
  reduction = "UMAP",
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunSlingshot-2.svg", width = 9, height = 3)


thisutils::log_message("=== Section: Monocle3 ===")

if (interactive()) {
  pancreas_sub <- RunMonocle3(
    pancreas_sub,
    group.by = "SubCellType"
  )
  trajectory <- pancreas_sub@tools$Monocle3$trajectory
  milestones <- pancreas_sub@tools$Monocle3$milestones
  p <- CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_partitions",
    reduction = "UMAP",
    label = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  ) +
    trajectory +
    milestones +
    CellDimPlot(
      pancreas_sub,
      group.by = "Monocle3_clusters",
      reduction = "UMAP",
      label = TRUE,
      palette = "ChineseSet8",
      xlab = "UMAP_1",
      ylab = "UMAP_2"
    ) +
    trajectory +
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = "UMAP",
      label = TRUE,
      palette = "ChineseSet8",
      xlab = "UMAP_1",
      ylab = "UMAP_2"
    ) +
    trajectory +
    FeatureDimPlot(
      pancreas_sub,
      features = "Monocle3_Pseudotime",
      reduction = "UMAP",
      xlab = "UMAP_1",
      ylab = "UMAP_2"
    ) +
    trajectory
  save_fig(p, "RunMonocle3.svg", width = 10, height = 6)
  thisutils::log_message(
    "  RunMonocle3-3.png — produced by RunMonocle3() internally."
  )
}

thisutils::log_message("=== Section: Dynamic features ===")

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200
)
pancreas_sub <- AnnotateFeatures(
  pancreas_sub,
  species = "Mus_musculus",
  db = c("TF", "CSPA")
)

ht <- DynamicHeatmap(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  use_fitted = TRUE,
  use_raster = TRUE,
  raster_by_magick = FALSE,
  ht_params = ht_raster_params,
  n_split = 6,
  reverse_ht = "Lineage1",
  species = "Mus_musculus",
  db = "GO_BP",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE,
  exp_legend_title = "Z-score",
  cell_annotation_palette = "ChineseSet8",
  heatmap_palette = "viridis",
  cell_annotation = "SubCellType",
  separate_annotation = list(
    "SubCellType",
    c("Nnat", "Irx1")
  ),
  separate_annotation_palette = c("ChineseSet8", "Set1"),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"),
    c("forestgreen")
  ),
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  height = 5,
  width = 2
)
if (just_margin) {
  p <- ht$plot &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
} else {
  p <- ht$plot
}
save_fig(
  p,
  "DynamicHeatmap-1.svg",
  width = 17.9,
  height = 6.8
)

p <- DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage1", "Lineage2"),
  group.by = "SubCellType",
  features = c(
    "Plk1",
    "Hes1",
    "Neurod2",
    "Ghrl",
    "Gcg",
    "Ins2"
  )
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "DynamicPlot-1.svg", width = 9.5, height = 5)

p <- FeatureStatPlot(
  pancreas_sub,
  group.by = "SubCellType",
  bg.by = "CellType",
  stat.by = c("Sox9", "Neurod2", "Isl1", "Rbp4"),
  add_box = TRUE,
  palette = "ChineseSet8",
  bg_palette = "ChineseSet8",
  comparisons = list(
    c("Ductal", "Ngn3 low EP"),
    c("Ngn3 high EP", "Pre-endocrine"),
    c("Alpha", "Beta")
  )
) +
  patchwork::plot_layout(
    guides = "collect"
  )
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "FeatureStatPlot-1.svg", width = 8, height = 5)

thisutils::log_message("=== Section: Differential expression analysis ===")

pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType",
  fc.threshold = 1,
  only.pos = FALSE
)

p <- DEtestPlot(
    pancreas_sub,
    group.by = "CellType",
    plot_type = "volcano",
  group_palette = "ChineseSet8",
  label.size = 2
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
p <- thisplot::rasterise_plot(p, dpi = raster_dpi)
save_fig(p, "RunDEtest-1.svg", width = 9, height = 4.5)

p <- DEtestPlot(
    pancreas_sub,
    group.by = "CellType",
    plot_type = "manhattan",
  group_palette = "ChineseSet8",
  pt.size = 1,
  label.size = 2.5
) +
  ggplot2::theme(aspect.ratio = 1 / 2)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
p <- thisplot::rasterise_plot(p, dpi = raster_dpi)
save_fig(p, "RunDEtest-2.svg", width = 9, height = 2.5)

p <- DEtestPlot(
    pancreas_sub,
    group.by = "CellType",
    plot_type = "ring",
  group_palette = "ChineseSet8",
  pt.size = 1,
  label.size = 2
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
p <- thisplot::rasterise_plot(p, dpi = raster_dpi)
save_fig(p, "RunDEtest-3.svg", width = 9, height = 3)

DEGs <- pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]

ht <- FeatureHeatmap(
  pancreas_sub,
  group.by = "CellType",
  features = DEGs$gene,
  feature_split = DEGs$group1,
  use_raster = TRUE,
  raster_by_magick = FALSE,
  ht_params = ht_raster_params,
  group_palette = "ChineseSet8",
  exp_legend_title = "Z-score",
  species = "Mus_musculus",
  db = c("GO_BP", "KEGG", "WikiPathway"),
  anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"),
    c("forestgreen")
  ),
  height = 5,
  width = 3
)
if (just_margin) {
  p <- ht$plot &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
} else {
  p <- ht$plot
}
save_fig(
  p,
  "FeatureHeatmap-1.svg",
  width = 15.2,
  height = 5.8
)


thisutils::log_message(
  "=== Section: Enrichment analysis (over-representation) ==="
)

pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group.by = "CellType",
  db = "GO_BP",
  species = "Mus_musculus",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "bar"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "RunEnrichment-1.svg",
  width = 9,
  height = 2.5
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "wordcloud"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "RunEnrichment-2.svg",
  width = 9,
  height = 3
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = c("Ductal", "Endocrine"),
  plot_type = "wordcloud",
  word_type = "feature"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "RunEnrichment-3.svg",
  width = 9,
  height = 3
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = "Ngn3-low-EP",
  plot_type = "network"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "RunEnrichment-4.svg",
  width = 9,
  height = 3.5
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = "Ductal",
  plot_type = "enrichmap"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "Enrichment_enrichmap-1.svg",
  width = 9,
  height = 6
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  plot_type = "comparison",
  topTerm = 3
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "Enrichment_comparison-1.svg",
  width = 10,
  height = 4
)

p <- EnrichmentPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = "Ductal",
  plot_type = "lollipop"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(
  p,
  "Enrichment_lollipop-1.svg",
  width = 8,
  height = 2.8
)


thisutils::log_message("=== Section: Enrichment analysis (GSEA) ===")

pancreas_sub <- RunGSEA(
  pancreas_sub,
  group.by = "CellType",
  db = "GO_BP",
  species = "Mus_musculus",
  DE_threshold = "p_val_adj < 0.05"
)

p <- GSEAPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = "Endocrine",
  id_use = "GO:0007186"
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "RunGSEA-1.svg", width = 8.5, height = 4)

p <- GSEAPlot(
  pancreas_sub,
  group.by = "CellType",
  group_use = "Endocrine",
  plot_type = "bar",
  direction = "both",
  topTerm = 10
)
if (just_margin) {
  p <- p &
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}
save_fig(p, "GSEA_bar-1.svg", width = 8.5, height = 3)


thisutils::log_message("=== Section: SCExplorer ===")
thisutils::log_message(
  "  [NOTE] SCExplorer-1.png is a screenshot of the interactive Shiny app.
        Run the following interactively to launch the app:
        PrepareSCExplorer(
        list(mouse_pancreas = pancreas_sub, human_pancreas = panc8_sub),
          base_dir = './SCExplorer'
        )
        app <- RunSCExplorer(base_dir = './SCExplorer')
        shiny::runApp(app)"
)


thisutils::log_message("=== Section: README composite examples ===")

thisutils::log_message("=== Section: CellDimPlot examples ===")

p_ex1_1 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE,
  legend.position = "none",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_2 <- CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  add_density = TRUE,
  label = TRUE,
  theme_use = "theme_blank",
  legend.position = "none",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_3 <- CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  stat.by = "Phase",
  reduction = "UMAP",
  stat_plot_type = "pie",
  stat_plot_size = 0.14,
  cells.highlight = TRUE,
  theme_use = "theme_blank",
  legend.position = "none",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_4 <- CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  hex = TRUE,
  hex.bins = 20,
  legend.position = "none",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_5 <- CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "UMAP",
  graph = "Standardpca_SNN",
  edge_color = "grey75",
  legend.position = "none",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_6 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_span = 0.1,
  legend.position = "none",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_7 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  lineages = paste0("Lineage", 1:2),
  lineages_whiskers = TRUE,
  legend.position = "none",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_8 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@misc$paga,
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE,
  legend.position = "none",
  theme_use = "theme_blank",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_9 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 3,
  pt.alpha = 0.3,
  velocity = "stochastic",
  velocity_plot_type = "grid",
  theme_use = "theme_blank",
  legend.position = "none",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1_10 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  pt.size = 3,
  pt.alpha = 0.3,
  velocity = "stochastic",
  velocity_plot_type = "stream",
  velocity_density = 2,
  velocity_smooth = 1,
  streamline_n = 20,
  theme_use = "theme_blank",
  legend.position = "none",
  palette = "ChineseSet8",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p_ex1 <- patchwork::wrap_plots(
  p_ex1_1,
  p_ex1_2,
  p_ex1_3,
  p_ex1_4,
  p_ex1_5,
  p_ex1_6,
  p_ex1_7,
  p_ex1_8,
  p_ex1_9,
  p_ex1_10,
  ncol = 5,
  byrow = TRUE
)
p_ex1 <- thisplot::rasterise_plot(p_ex1, dpi = raster_dpi)
save_fig(
  p_ex1,
  "Example-1.svg",
  width = 15,
  height = 6,
  dpi = examples_dpi
)

thisutils::log_message("=== Section: CellStatPlot examples ===")

p_ex2_1 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "SubCellType",
  bg.by = "CellType",
  palette = "Set1",
  stat_type = "count",
  position = "dodge"
)

p_ex2_2 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  position = "stack",
  palette = "Set2",
  plot_type = "bar"
)

p_ex2_3 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  palette = "Spectral",
  position = "dodge",
  plot_type = "rose"
)

p_ex2_4 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  position = "dodge",
  palette = "RdYlBu",
  plot_type = "ring"
)

p_ex2_12 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "percent",
  position = "stack",
  palette = "Reds",
  plot_type = "bar"
)

p_ex2_13 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "percent",
  position = "stack",
  palette = "Blues",
  plot_type = "area"
)

p_ex2_5 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "percent",
  palette = "Dark2",
  plot_type = "dot"
)

p_ex2_6 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "count",
  position = "stack",
  palette = "ChineseSet8",
  plot_type = "trend"
)

p_ex2_6_2 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  group.by = "CellType",
  stat_type = "percent",
  position = "stack",
  palette = "Set3",
  plot_type = "trend"
)

p_ex2_7 <- CellStatPlot(
  pancreas_sub,
  stat.by = "Phase",
  plot_type = "pie",
  palette = "simspec"
)

p_ex2_8 <- CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  plot_type = "sankey",
  palette = "Geyser"
)

p_ex2_9 <- CellStatPlot(
  pancreas_sub,
  stat.by = c("CellType", "Phase"),
  legend.position = "none",
  plot_type = "chord"
)

p_ex2_9 <- if (ggplot2::is.ggplot(p_ex2_9) || inherits(p_ex2_9, "grob")) {
  p_ex2_9
} else {
  temp_png <- tempfile(fileext = ".png")
  on.exit(unlink(temp_png), add = TRUE)
  grDevices::png(filename = temp_png, width = 1800, height = 1400, res = 200)
  CellStatPlot(
    pancreas_sub,
    stat.by = c("CellType", "Phase"),
    legend.position = "none",
    plot_type = "chord"
  )
  grDevices::dev.off()
  if (!requireNamespace("png", quietly = TRUE)) {
    log_message(
      "Package 'png' is required to compose chord plot in patchwork",
      message_type = "error"
    )
  }
  if (!file.exists(temp_png) || file.info(temp_png)$size <= 0) {
    log_message(
      "Failed to render chord plot for patchwork composition",
      message_type = "error"
    )
  }
  captured <- png::readPNG(temp_png)
  patchwork::wrap_elements(
    full = grid::rasterGrob(captured, interpolate = TRUE)
  )
}

pancreas_sub$G2M <- pancreas_sub$Phase == "G2M"
pancreas_sub$S <- pancreas_sub$Phase == "S"
pancreas_sub$G1 <- pancreas_sub$Phase == "G1"

p_ex2_10 <- CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "S",
    "G2M",
    "G1"
  ),
  plot_type = "venn",
  stat_level = "TRUE",
  palette = "ChineseSet8"
)

p_ex2_11 <- CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "S",
    "G2M",
    "G1"
  ),
  plot_type = "upset",
  stat_level = "TRUE"
)

p_ex2_top <- patchwork::wrap_plots(
  p_ex2_1,
  p_ex2_2,
  p_ex2_12,
  p_ex2_13,
  p_ex2_6,
  p_ex2_6_2,
  p_ex2_5,
  p_ex2_3,
  p_ex2_4,
  p_ex2_7,
  ncol = 5,
  byrow = TRUE
)

p_ex2_bottom <- patchwork::wrap_plots(
  p_ex2_8,
  p_ex2_10,
  p_ex2_11,
  p_ex2_9,
  ncol = 4,
  widths = c(0.3, 0.3, 0.2, 0.2),
  byrow = TRUE
)

p_ex2 <- patchwork::wrap_plots(
  p_ex2_top,
  p_ex2_bottom,
  ncol = 1,
  heights = c(0.7, 0.3),
  byrow = TRUE
)
save_fig(
  p_ex2,
  "Example-2.svg",
  width = 19,
  height = 11,
  dpi = examples_dpi
)


thisutils::log_message("=== Section: FeatureStatPlot examples ===")
p_ex3_1 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Fev",
  group.by = "SubCellType",
  palette = "ChineseSet8",
  plot_type = "violin"
)

p_ex3_2 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Fev",
  group.by = "SubCellType",
  palette = "ChineseSet8",
  plot_type = "bar"
)

p_ex3_3 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Rbp4",
  group.by = "SubCellType",
  palette = "ChineseSet8",
  plot_type = "box"
)

p_ex3_4 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Rbp4",
  group.by = "SubCellType",
  palette = "ChineseSet8",
  plot_type = "dot"
)

p_ex3_5 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Pyy",
  group.by = "SubCellType",
  palette = "ChineseSet8",
  add_box = TRUE,
  comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta"))
)

p_ex3_7 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Fev",
  group.by = "SubCellType",
  fill.by = "expression",
  palette = "ChineseSet8",
  bg_palette = "ChineseSet8",
  bg.by = "CellType",
  same.y.lims = TRUE
)

p_ex3_8 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "Rbp4",
  group.by = "SubCellType",
  palette = "ChineseSet8",
  plot_type = "col"
)

p_ex3_9 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = c(
    "Sox9",
    "Anxa2",
    "Bicc1",
    "Neurog3",
    "Hes6",
    "Fev",
    "Neurod1",
    "Rbp4",
    "Pyy",
    "Ins1",
    "Gcg",
    "Sst",
    "Ghrl"
  ),
  group.by = "SubCellType",
  plot.by = "feature",
  palette = "ChineseSet8",
  stack = TRUE
)

p_ex3_10 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = "G2M_score",
  group.by = "SubCellType",
  split.by = "Phase",
  plot_type = "violin",
  palette = "ChineseSet8",
  comparisons = TRUE
)

p_ex3_11 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = c(
    "Sox9",
    "Anxa2",
    "Bicc1",
    "Neurog3",
    "Hes6",
    "Fev",
    "Neurod1",
    "Rbp4",
    "Pyy",
    "Ins1",
    "Gcg",
    "Sst",
    "Ghrl"
  ),
  fill.by = "feature",
  plot_type = "box",
  group.by = "SubCellType",
  bg.by = "CellType",
  # flip = TRUE,
  palette = "ChineseSet8",
  stack = TRUE
)

p_ex3_12 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = c(
    "Sox9",
    "Anxa2",
    "Bicc1",
    "Neurog3",
    "Hes6",
    "Fev",
    "Neurod1",
    "Rbp4",
    "Pyy",
    "Ins1",
    "Gcg",
    "Sst",
    "Ghrl"
  ),
  group.by = "SubCellType",
  bg.by = "CellType",
  palette = "ChineseSet8",
  stack = TRUE
)

p_ex3_13 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = c("G2M_score"),
  group.by = "SubCellType",
  palette = "ChineseSet8",
  add_point = TRUE
)

p_ex3_14 <- FeatureStatPlot(
  pancreas_sub,
  stat.by = c("Fev"),
  group.by = "SubCellType",
  palette = "ChineseSet8",
  add_trend = TRUE
)

p_ex3_p2_left <- patchwork::wrap_plots(
  p_ex3_13,
  p_ex3_14,
  ncol = 1,
  guides = "collect",
  byrow = TRUE
)

p_ex3_p1 <- patchwork::wrap_plots(
  p_ex3_1,
  p_ex3_2,
  p_ex3_7,
  p_ex3_3,
  p_ex3_4,
  p_ex3_8,
  p_ex3_5,
  p_ex3_10,
  ncol = 4,
  byrow = TRUE
)
p_ex3_p2 <- patchwork::wrap_plots(
  p_ex3_p2_left,
  p_ex3_9,
  p_ex3_11,
  p_ex3_12,
  ncol = 4,
  widths = c(0.2, 0.4, 0.3, 0.2),
  byrow = TRUE
)
p_ex3 <- patchwork::wrap_plots(
  p_ex3_p1,
  p_ex3_p2,
  ncol = 1,
  heights = c(0.55, 0.45),
  byrow = TRUE
)
save_fig(
  p_ex3,
  "Example-3.svg",
  width = 18,
  height = 12,
  dpi = examples_dpi
)


thisutils::log_message("=== Section: GroupHeatmap examples ===")
features_ex4 <- c(
  "Sox9",
  "Anxa2",
  "Bicc1",
  "Neurog3",
  "Hes6",
  "Fev",
  "Neurod1",
  "Rbp4",
  "Pyy",
  "Ins1",
  "Gcg",
  "Sst",
  "Ghrl"
)
ht_ex4_1 <- GroupHeatmap(
  pancreas_sub,
  ht_params = ht_raster_params,
  features = features_ex4,
  nlabel = 0,
  exp_legend_title = "Z-score",
  width = 1.6,
  height = 2.8,
  group_palette = "ChineseSet8",
  heatmap_palette = "viridis",
  show_row_names = TRUE,
  group.by = c("CellType", "SubCellType")
)

pancreas_sub <- RunDEtest(
  pancreas_sub,
  group.by = "CellType"
)
de_filter <- dplyr::filter(
  pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
  p_val_adj < 0.05 & avg_log2FC > 1
)

ht_ex4_2 <- GroupHeatmap(
  pancreas_sub,
  ht_params = ht_raster_params,
  nlabel = 0,
  exp_legend_title = "Z-score",
  features = de_filter$gene,
  group.by = "SubCellType",
  split.by = "Phase",
  # n_split = 2,
  # anno_terms = TRUE,
  # terms_width = 0.3,
  # cores = 2,
  width = 3.5,
  height = 2.6,
  group_palette = "ChineseSet8",
  cell_split_palette = "Dark2"
)

ht_ex4_3 <- GroupHeatmap(
  pancreas_sub,
  ht_params = ht_raster_params,
  exp_legend_title = "Z-score",
  features = features_ex4,
  group.by = "SubCellType",
  add_violin = TRUE,
  cluster_rows = TRUE,
  fill.by = "expression",
  group_palette = "ChineseSet8",
  fill_palette = "RdYlBu",
  nlabel = 0,
  show_row_names = TRUE,
  width = 3,
  height = 2.8
)

ht_ex4_4 <- GroupHeatmap(
  pancreas_sub,
  ht_params = ht_raster_params,
  exp_legend_title = "Z-score",
  features = features_ex4,
  group.by = "CellType",
  split.by = "Phase",
  n_split = 4,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  add_dot = TRUE,
  dot_size = grid::unit(6, "mm"),
  add_reticle = TRUE,
  group_palette = "ChineseSet8",
  heatmap_palette = "viridis",
  nlabel = 0,
  show_row_names = TRUE,
  width = 2.8,
  height = 2.6
)

de_top <- de_filter |>
  dplyr::group_by(gene) |>
  dplyr::top_n(1, avg_log2FC) |>
  dplyr::group_by(group1) |>
  dplyr::top_n(3, avg_log2FC)
ht_ex4_5 <- GroupHeatmap(
  pancreas_sub,
  ht_params = ht_raster_params,
  exp_legend_title = "Z-score",
  features = de_top$gene,
  feature_split = de_top$group1,
  group.by = "CellType",
  heatmap_palette = "YlOrRd",
  cell_annotation = c(
    "Phase",
    "G2M_score",
    "Neurod2"
  ),
  cell_annotation_palette = c(
    "Dark2",
    "ChineseSet8",
    "ChineseSet8"
  ),
  cell_annotation_params = list(
    width = grid::unit(5, "mm")
  ),
  row_title_rot = 90,
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"),
    c("forestgreen")
  ),
  add_dot = TRUE,
  dot_size = grid::unit(6, "mm"),
  add_bg = TRUE,
  nlabel = 0,
  flip = TRUE,
  width = 2.7,
  height = 1.8,
  show_row_names = FALSE,
  show_column_names = FALSE
)

ht_ex4_6 <- GroupHeatmap(
  pancreas_sub,
  ht_params = ht_raster_params,
  exp_legend_title = "Z-score",
  features = features_ex4,
  group.by = c("CellType", "SubCellType"),
  add_dot = TRUE,
  dot_size = grid::unit(6, "mm"),
  group_palette = "ChineseSet8",
  heatmap_palette = "plasma",
  add_reticle = TRUE,
  show_row_names = TRUE,
  nlabel = 0,
  width = 1.6,
  height = 3.2
)

p_ex4 <- patchwork::wrap_plots(
  ht_ex4_2$plot,
  ht_ex4_1$plot,
  ht_ex4_3$plot,
  ht_ex4_4$plot,
  ht_ex4_6$plot,
  ht_ex4_5$plot,
  ncol = 3,
  byrow = TRUE
)
save_fig(
  p_ex4,
  "Example-4.svg",
  width = 20,
  height = 7.6,
  dpi = examples_dpi
)


thisutils::log_message("=== Section: Raster-only figures ===")
for (figure in external_figures) {
  save_dark_png(figure)
}

thisutils::log_message("========================================")
thisutils::log_message("Figure reproduction complete!")
thisutils::log_message("========================================")
thisutils::log_message("Output directory: {.file {fig_dir}}")

variants <- c(".svg", "-dark.svg")
expected_figures <- c(
  "EDA-1",
  "EDA-2",
  "EDA-3",
  "EDA-4",
  "EDA-5",
  "RunCellQC-1",
  "RunCellQC-2",
  "RunCellQC-3",
  "integration_scop-1",
  "RunKNNMap-1",
  "RunKNNPredict-bulk-1",
  "RunKNNPredict-scrna-1",
  "RunKNNPredict-scrna-3",
  "RunCytoTRACE",
  "RunSCVELO-1",
  "RunSCVELO-2",
  "RunPAGA-1",
  "RunSlingshot-1",
  "RunSlingshot-2",
  "RunMonocle3",
  "DynamicHeatmap-1",
  "DynamicPlot-1",
  "FeatureStatPlot-1",
  "RunDEtest-1",
  "RunDEtest-2",
  "RunDEtest-3",
  "FeatureHeatmap-1",
  "RunEnrichment-1",
  "RunEnrichment-2",
  "RunEnrichment-3",
  "RunEnrichment-4",
  "Enrichment_enrichmap-1",
  "Enrichment_comparison-1",
  "Enrichment_lollipop-1",
  "RunGSEA-1",
  "GSEA_bar-1",
  "SCExplorer-1",
  "Example-1",
  "Example-2",
  "Example-3",
  "Example-4"
)

expected_files <- unlist(lapply(expected_figures, function(figure) {
  if (figure %in% external_figures) {
    paste0(figure, ".png")
  } else {
    paste0(figure, variants)
  }
}))
raster_files <- paste0(external_figures, "-dark.png")
existing <- file.exists(file.path(fig_dir, c(expected_files, raster_files)))
status <- function(present) {
  ifelse(is.na(present), "-", ifelse(present, "OK", "MISSING"))
}

thisutils::log_message(sprintf(
  "  %-35s %-6s %-9s",
  "Figure",
  ".svg",
  "-dark.svg"
))
thisutils::log_message(sprintf(
  "  %-35s %-6s %-9s",
  "------",
  "----",
  "---------"
))
for (figure in setdiff(expected_figures, external_figures)) {
  present <- file.exists(file.path(fig_dir, paste0(figure, variants)))
  thisutils::log_message(sprintf(
    "  %-35s %-6s %-9s",
    figure,
    status(present[1]),
    status(present[2])
  ))
}
thisutils::log_message(
  "  Raster-only figures (dark copy is {.file -dark.png}):"
)
thisutils::log_message(sprintf(
  "  %-35s %-6s %-9s",
  "Figure",
  ".png",
  "-dark.png"
))
thisutils::log_message(sprintf(
  "  %-35s %-6s %-9s",
  "------",
  "----",
  "---------"
))
for (i in seq_along(external_figures)) {
  thisutils::log_message(sprintf(
    "  %-35s %-6s %-9s",
    external_figures[i],
    status(file.exists(file.path(
      fig_dir,
      paste0(external_figures[i], ".png")
    ))),
    status(file.exists(file.path(fig_dir, raster_files[i])))
  ))
}
thisutils::log_message(
  "Total: {.val {sum(existing)}} / {.val {length(existing)}} figure files present"
)
thisutils::log_message(
  "Figures are written straight into {.file {fig_dir}}; commit and push the figures repository to update the README."
)
