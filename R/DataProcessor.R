# For data processing
library(devtools)
library(dplyr)
library(Matrix)
library(Seurat)

# For data visualizing
library(concaveman)
library(ggforce)
library(ggplot2)
library(ggsci)
library(graphics)
library(patchwork)

# Load helper functions in other files
source('R/Utility.R')
source('R/PlotStyle.R')

# Create a results folder
CreateResultsFolder('Results')

csf_data <- read.table(
  gzfile('Data/csf_data.csv.gz'), 
  header = T, 
  row.names = 1, 
  sep = '\t'
)
csf_data <- log(csf_data + 1)

# Load and process single-cell dataset
cell_types <- unlist(lapply(colnames(csf_data), ExtractField, 1))
table(cell_types)

# Create seurat object. Ensure all genes have at least 3 cells
# and every cell has at least 200 features
csf <- Seurat::CreateSeuratObject(
  counts = csf_data,
  min.cells = 3,
  min.features = 200,
  project = 'CSF_2020')

# AddMetaData adds columns to object@meta.data:
sc_cell_info <- utils::read.table('Data/single_cell_metadata.txt', header = T, row.names=1)
csf <- Seurat::AddMetaData(
  object = csf,
  metadata = sc_cell_info,
  col.name = c('Twin', 'Case', 'Sample', 'index.sort', 'Clones')
)

cell_info_matrix = Matrix::as.matrix(sc_cell_info)
# Get the subset of `sc_cell_info` whose `Sample` is 'CSF'
cell_info_csf = subset(sc_cell_info, sc_cell_info[, 'Sample'] == 'CSF')
# Replace all 'CD4Tcell' in `index.sort` with 'CD4'
cell_info_csf$index.sort[cell_info_csf$index.sort == 'CD4Tcell'] <- 'CD4'
# Replace all 'CD8Tcell' in `index.sort` with 'CD8'
cell_info_csf$index.sort[cell_info_csf$index.sort == 'CD8Tcell'] <- 'CD8'

# A. Process data whose `index.sort` is 'CD4'
cell_info_cd4 <- subset(cell_info_csf, cell_info_csf[, 'index.sort'] == 'CD4')

# A-1. Get HD only subset
cell_info_cd4_hd <- subset(cell_info_cd4, cell_info_cd4[, 'Case'] == 'HD')

distribution_cd4_hd <-
  ggplot(cell_info_cd4_hd, aes(x = Twin, fill = Clones)) +
  geom_bar(stat = 'count', position = 'stack', width = 0.7) +
  geom_text(aes(label = ..count..), stat = 'count', size = 3.5, position = position_stack(vjust = .5)) +
  StandardPlotTheme +
  scale_fill_manual(values = StandardFills) +
  labs(title = "Distribution of CD4 in HD") +
  xlab(label = "Twin") +
  ylab(label = "Count")

ggsave(
  file = 'Results/dis_cd4_with_hd.pdf',
  plot = distribution_cd4_hd,
  width = 5,
  height = 4
)

proportion_cd4_hd <-
  ggplot(cell_info_cd4_hd, aes(x = Twin, fill = Clones)) +
  geom_bar(stat = 'count', position = 'fill', width = 0.7) +
  geom_text(
    # Use `ase(label = ..count..)` to label the total count
    aes(label = scales::percent(..count../sum(..count..))),
    stat = 'count',
    size = 3.5,
    position = position_fill(vjust = .5)
  ) +
  StandardPlotTheme +
  scale_fill_manual(values = StandardFills) +
  labs(title = "Proportion of CD4 in HD") +
  xlab(label = "Twin") +
  ylab(label = "Count")

ggsave(
  file = 'Results/pro_cd4_with_hd.pdf',
  plot = proportion_cd4_hd,
  width = 5,
  height = 4
)

# A-2. Get MS only subset

cell_info_cd4_ms <- subset(cell_info_cd4, cell_info_cd4[, 'Case'] == 'MS')

distribution_cd4_ms <-
  ggplot(cell_info_cd4_ms, aes(x = Twin, fill = Clones)) +
  geom_bar(stat = 'count', position = 'stack', width = 0.7) +
  geom_text(aes(label = ..count..), stat = 'count', size = 3.5, position = position_stack(vjust = .5)) +
  StandardPlotTheme +
  scale_fill_manual(values = StandardFills) +
  labs(title = "Distribution of CD4 in MS") +
  xlab(label = "Twin") +
  ylab(label = "Count")

ggsave(
  file = 'Results/dis_cd4_with_ms.pdf',
  plot = distribution_cd4_ms,
  width = 5,
  height = 4
)

proportion_cd4_ms <-
  ggplot(cell_info_cd4_ms, aes(x = Twin, fill = Clones)) +
  geom_bar(stat = 'count', position = 'fill', width = 0.7) +
  geom_text(
    # Use `ase(label = ..count..)` to label the total count
    aes(label = scales::percent(..count../sum(..count..))),
    stat = 'count',
    size = 2,
    position = position_fill(vjust = .5)
  ) +
  StandardPlotTheme +
  scale_fill_manual(values = StandardFills) +
  labs(title = "Proportion of CD4 in MS") +
  xlab(label = "Twin") +
  ylab(label = "Count")

ggsave(
  file = 'Results/pro_cd4_with_ms.pdf',
  plot = proportion_cd4_ms,
  width = 5,
  height = 4
)

# B. Process data whose `index.sort` is 'CD8'
cell_info_cd8 <- subset(cell_info_csf, cell_info_csf[, 'index.sort'] == 'CD8')

# QC and selecting cells for further analysis.
# We calculate the percentage of mitochondrial genes here
# and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and non-log-normalized counts.
# The % of counts mapping to MT-genes is a common scRNAseq QC metric.

mito.genes <- grep(
  pattern = "^MT-",
  x = rownames(x = GetAssayData(csf)),
  value = TRUE)

percent.mito <-
  Matrix::colSums(GetAssayData(csf, slot = 'counts')[mito.genes, ]) / Matrix::colSums(GetAssayData(csf, slot = 'counts'))

csf <- AddMetaData(object = csf, metadata = percent.mito, col.name = 'percent.mito')

VlnPlot(
  object = csf,
  features = c('nFeature_RNA', 'percent.mito'),
  ncol = 2,
  slot = 'counts',
  flip = TRUE
)

csf <- subset(csf, subset = percent.mito > -Inf & percent.mito < 0.05)
csf <- subset(csf, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
#csf <- subset(csf, subset = CD4 > 0 | CD8A > 0)
csf <- subset(
  csf,
  subset =
    index.sort == 'CD8' |
    index.sort == 'CD4' |
    index.sort == 'CD8Tcell' |
    index.sort == 'CD4Tcell'
)
csf <- subset(csf, subset = Sample == 'CSF')
csf <- subset(csf, subset = Case == 'MS' | Case == 'HD')

# Normalizing the data
# Further normalization is performed since the dataset used to create
# the Seurat object is a merge of normalized datasets into one file.
# That means that, in the final merged file, the data from different dataset
# are not normalized to the same sequencing depth.
# Therefore, a further normalization is required.

csf <- NormalizeData(object = csf, normalization.method = 'LogNormalize', scale.factor = 10000)

# Identification of highly variable features

csf <- FindVariableFeatures(
  csf,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  mean.cutoff = c(0.0125, 3),
  dispersion.cutoff = c(0.5,Inf)
)

# Number of highly variable features
length(x = VariableFeatures(csf))

# Scaling the data and removing unwanted sources of variation
csf <- ScaleData(csf, vars.to.regress = 'percent.mito')

# Print the top 10 highly variable features
head(VariableFeatures(csf), 10)

# Draw variable features
VariableFeaturePlot(csf)

# Perform linear dimensional reduction
csf <- RunPCA(csf, features = VariableFeatures(csf), verbose = TRUE, ndims.print = 1:5, nfeatures.print = 5)

# Visulize the 1st and 2nd principle components
#VizDimLoadings(csf, dims = 1:2, reduction = "pca")

# Visualize PCA plot
#DimPlot(csf, reduction = "pca", label = FALSE, dims = c(1, 2))

# ProjectDim scores each gene in the dataset (including genes not included in the PCA)
# based on their correlation with the calculated components

# Project Dimensional Reduction Onto Full CSF dataset
# And print the 1st and 2nd principle components
csf <- ProjectDim(csf, reduction = "pca", dims.print = 1:5, nfeatures.print = 5, verbose = TRUE)

# Heatmaps based on the PCA
DimHeatmap(csf, dims = 1:2, cells = 500, balanced = TRUE)

hdOnlyCSF <- subset(csf, subset = Case == 'HD')
DimHeatmap(hdOnlyCSF, dims = 1:2, cells = 500, balanced = TRUE)

# Determine statistically significant principal components
csf <- JackStraw(csf, num.replicate = 100)
csf <- ScoreJackStraw(csf, dims = 1:5)
JackStrawPlot(csf, dims = 1:5, reduction = 'pca')
# Quickly Pick Relevant Dimensions
ElbowPlot(csf)

# Cluster the cells
csf <- FindNeighbors(csf, reduction = 'pca', dims = 1:5, graph.name = 'snn_graph')
csf <- FindClusters(csf, graph.name = 'snn_graph', resolution = 1, verbose = TRUE)

# Run Non-linear dimensional reduction (tSNE)
csf <- RunTSNE(
  csf,
  dims = 1:10,
  do.fast = TRUE,
  perplexity = 10,
  reduction.key = 'tSNE_',
  check_duplicates = FALSE
)

dims_cd4_cd8 <- DimPlot(
  csf,
  reduction = "tsne",
  group.by = 'index.sort',
  cols = c("#00c3e3", "#00c3e3", "#f70d1a"),
  #order =  c("CD4","CD4Tcell","CD8Tcell"),
  pt.size = 1,
  label.size = 3,
  label = FALSE
) +
  StandardPlotTheme +
  labs(title = "Dim: CD4 X CD8") +
  xlab(label = "tSNE1") +
  ylab(label = "tSNE2")

ggsave(
  file = 'Results/dim_cd4_with_cd8.pdf',
  plot = dims_cd4_cd8,
  width = 5,
  height = 4
)

dims_cd4_cd8_split <- DimPlot(
  csf,
  reduction = "tsne",
  group.by = 'index.sort',
  split.by = 'Case',
  cols = c("#00c3e3", "#00c3e3", "#f70d1a"),
  #order =  c("CD4","CD4Tcell","CD8Tcell"),
  pt.size = 1,
  label.size = 3,
  label = FALSE
) +
  StandardPlotTheme +
  labs(title = "Dim CD4 X CD8 - Split by HD & MS") +
  xlab(label = "tSNE1") +
  ylab(label = "tSNE2")

ggsave(
  file = 'Results/dim_cd4_with_cd8_split_by_cases.pdf',
  plot = dims_cd4_cd8_split,
  width = 8,
  height = 4
)

DimPlot(
  csf,
  reduction = "tsne",
  group.by = 'Case',
  label.size = 3,
  label = TRUE
)

plots <- FeaturePlot(
  csf,
  features = c('CD4', 'CD8A'),
  cols = c("#f0f0f0", "#00c3e3", "#f70d1a"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  reduction = 'tsne',
  blend = TRUE,
  blend.threshold = 0.5
)

features_cd4_cd8a <- plots[[3]] +
  StandardPlotTheme +
  labs(title = 'Feature: CD4 X CD8A') +
  xlab(label = 'tSNE1') +
  ylab(label = 'tSNE2')

color_blend_info <- plots[[4]] +
  StandardPlotTheme +
  labs(title = 'Color Blend: 0.5') +
  xlab(label = 'CD4') +
  ylab(label = 'CD8A')

ggsave(
  file = 'Results/feature_cd4_with_cd8a.pdf',
  plot = features_cd4_cd8a + color_blend_info,
  width = 8,
  height = 4
)
