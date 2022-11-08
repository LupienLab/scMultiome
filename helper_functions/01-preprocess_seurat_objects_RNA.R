# This script contains helper functions to preprocess seurat object for the RNA

library(Seurat)
library(Matrix)

# Adapted from Signac with modifications
# Written by: Shalini Bahl
# Contact: Shalini.Bahl@uhnresearch.ca

# VERSION: 1
# Last updated: SEPT 6, 2022

############## REFERENCE DATA #######################

# Note to self - transfer to a reference data directory here 
seqinfo <- readRDS('/cluster/projects/lupiengroup/Shalini/CommonData/hg38_seqinfo.rds')


############## FUNCTIONS ############################

#' Create a RNA counts matrix from the 10X multiome data
#'
#' @param Path to the cellranger-arc output directory
#'

test_generate_RNA_counts_matrix <- function(matrix_dir){
	inputdata.10x <- Read10X_h5(paste0(matrix_dir, "/outs/filtered_feature_bc_matrix.h5"))
	inputdata.10x <- inputdata.10x$`Gene Expression`
	return(inputdata.10x)
}

#' Add all the QC information to the RNA seurat object
#'
#' @param RNA seurat object
#'

add_RNA_QCs <- function(RNA_obj){
	### Add mitochrondrial information
	RNA_obj[["percent.mt"]] <- PercentageFeatureSet(RNA_obj, pattern = "^MT-")

	### Add the HB information
	HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
	HB_m <- match(HB.genes, rownames(RNA_obj@assays$RNA))	
	HB.genes <- rownames(RNA_obj@assays$RNA)[HB_m]
	HB.genes <- HB.genes[!is.na(HB.genes)]
	RNA_obj[["percent.HB"]] <- PercentageFeatureSet(RNA_obj, features=HB.genes)

	return(RNA_obj)
}

#' Compute the cell cycle scoring
#'
#' @param
#'

cell_cycling_qc <- function(RNA_obj){
	s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
        RNA_obj <- CellCycleScoring(RNA_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

	return(RNA_obj)
}

#' Filter RNA object after assigning optimal thresholds
#'
#' @param Seurat object for the RNA-seq
#' @param Mitochrondrial threshold for filtering
#' @param Red blood cells threshold for filtering
#' @param Number of genes threshold for filtering
#' @param Number of transcripts threshold for filtering


filter_RNA_obj <- function(RNA_obj, mito_thres, HB_thres, nFeature_thres, nCount_thres){
	RNA_obj <- subset(RNA_obj, subset = nFeature_RNA > nFeature_thres & nCount_RNA > nCount_thres & percent.mt < mito_thres & percent.HB < HB_thres)
	return(RNA_obj)
}

