# This script is to create unfiltered seurat objects and generate QC plots for manual assessment

library(tidyverse)
library(Seurat)
library(argparse)

source("./helper_functions/01-preprocess_seurat_objects_RNA.R")

#### Collect all the arguments ####################################################################

parser <- ArgumentParser()

parser$add_argument("-m", "--mat_dir", required=TRUE, help="Path to the cellranger-arc filtered_bc directory")
parser$add_argument("-s", "--samp", required=TRUE, help="New sample ID")
parser$add_argument("-col", "--sample_col", required=TRUE, help="HEX code for the sample colour")
parser$add_argument("-mito", "--mito_thresh", required=FALSE, help="Mitochondrial threshold for filtering")
parser$add_argument("-HB", "--HB_thresh", required=FALSE, help="Human blood cells threshold for filtering")
parser$add_argument("-gene", "--nFeature_thresh", required=FALSE, help="Gene threshold for filtering")
parser$add_argument("-count", "--nCount_thresh", required=FALSE, help="Transcript threshold for filtering")

args <- parser$parse_args()
print(args)

###################################################################################################

# Create the individual unfiltered seurat objects if they don't exist

unfiltered_seurat_obj <- paste0('./seurat_objects/unfiltered_seurat_obj/', args$samp, '_unfiltered_seurat_obj.rds')

if(file.exists(unfiltered_seurat_obj)){
RNA_obj <- readRDS(unfiltered_seurat_obj)
print(RNA_obj)
print(head(RNA_obj@meta.data))

} else {
## Generate the RNA counts matrix

RNA_counts <- test_generate_RNA_counts_matrix(args$mat_dir)

print(head(rownames(RNA_counts)))

## Create unfiltered seurat object
RNA_obj <- CreateSeuratObject(counts = RNA_counts, project = args$samp, min.cells = 1, min.features = 1) 
print(RNA_obj)

# QC stats
RNA_obj <- add_RNA_QCs(RNA_obj)

# Run normalization
RNA_obj <- NormalizeData(RNA_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Compute cell cycle scoring
RNA_obj <- cell_cycling_qc(RNA_obj)

print(head(RNA_obj@meta.data))

saveRDS(RNA_obj, paste0('./seurat_objects/unfiltered_seurat_obj/', args$samp, '_unfiltered_seurat_obj.rds'))
}

#### Generate relevant plots at the sample level - save the meta data so plots can be reformated later
qc_plots_dir <- paste0('./plots/QC/', args$samp)

Idents(RNA_obj) <- "orig.ident"

pdf(paste0(qc_plots_dir, "/", args$samp, "_feature_count_mito_pre_filtered.pdf"), width = 8)
VlnPlot(RNA_obj, features = c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.HB"), ncol = 4, pt.size = 0, col = args$sample_co) + NoLegend()
dev.off()

write.table(RNA_obj@meta.data, paste0(qc_plots_dir, "/", args$samp, "_qc_meta_data.tsv"), col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE)

#### Filter the individual seurat objects based on given parameters (needs to be consistent across all samples)
filtered_seurat_obj <- paste0('./seurat_objects/filtered_seurat_obj/', args$samp, '_filtered_seurat_obj.rds')

if(file.exists(filtered_seurat_obj)){
print("RNA seurat object already filtered!")} else{

if(args$mito_thres !="" & args$HB_thres != "" & args$nFeature_thres != "" & args$nCount_thres != ""){
print("Filtering seurat object and saving it")

RNA_obj <- filter_RNA_obj(RNA_obj, args$mito_thres, args$HB_thres, args$nFeature_thres, args$nCount_thres)
print(RNA_obj)

saveRDS(RNA_obj, paste0('./seurat_objects/filtered_seurat_obj/', args$samp, '_filtered_seurat_obj.rds'))

} else {print("Define the filtering parameters for analysis")}
}

sessionInfo()

