# This script is to generate some QC plots from merged samples for easy comparision

library(tidyverse)
library(Seurat)
library(argparse)

source("./helper_functions/merging_seurat_helper.R")

#### Collect all the arguments ####################################################################

parser <- ArgumentParser()

parser$add_argument("-f", "--meta_data", required=TRUE, help="The path to the .csv meta-data file")
parser$add_argument("-r", "--reads_cutoff", required=TRUE, help="Explicitly state what the reads cutoff is (added later)")
parser$add_argument("-gr", "--group", required=TRUE, help="Group name for samples to be merged")
parser$add_argument("-obj", "--obj_processed", required=TRUE, help="Should the 'filtered' or 'unfiltered' object be merged")

print(parser)

args <- parser$parse_args()
print(args)

###################################################################################################

print(args$reads_cutoff)

# Read in the meta-data containing custom colour information
meta_data <- read.table(args$meta_data, header = TRUE, sep = ',')
print(head(meta_data))

# Read in the merged seurat object to generate QC plots
merged_obj <- readRDS(paste0("./seurat_objects/", args$obj_processed, "_seurat_obj/merged_seurat_obj/", args$gr, "_", args$obj_processed, "_merged_seurat_obj.rds"))
print(merged_obj)
print(head(merged_obj@meta.data))

qc_plots_dir <- paste0('./plots/QC/', args$gr)

# QCs to plot: nFeature (number of genes), nCount (Number of transcripts)
transcripts_qc_plot <- generate_violin_plot(merged_obj, "nCount_RNA", args$reads_cutoff, "Number of transcripts")
#transcripts_qc_plot <- generate_violin_plot(merged_obj, "nCount_RNA", unique(meta_data$nCount_threshold), "Number of transcripts")
genes_qc_plot <- generate_violin_plot(merged_obj, "nFeature_RNA", unique(meta_data$nFeature_threshold), "Number of unique genes")
mito_qc_plot <- generate_violin_plot(merged_obj, "percent.mt", unique(meta_data$mito_threshold), "Percentage of reads mapping to mitochondrial genome")
HB_qc_plot <- generate_violin_plot(merged_obj, "percent.HB", unique(meta_data$HB_threshold), "Percentage of reads mapping to HB genes")

# Save the plots in a single PDF file (for now... TIFF format is preferred)
pdf(paste0(qc_plots_dir, "/", args$gr, "_", args$obj_processed, "_feature", unique(meta_data$nFeature_threshold), "_count", args$reads_cutoff, "_mito", unique(meta_data$mito_threshold), "_HB", unique(meta_data$HB_threshold), ".pdf"), width = 8)
#pdf(paste0(qc_plots_dir, "/", args$gr, "_", args$obj_processed, "_feature", unique(meta_data$nFeature_threshold), "_count", unique(meta_data$nCount_threshold), "_mito", unique(meta_data$mito_threshold), "_HB", unique(meta_data$HB_threshold), ".pdf"), width = 8)
print(transcripts_qc_plot)
print(genes_qc_plot)
print(mito_qc_plot)
print(HB_qc_plot)
dev.off()

sessionInfo()

