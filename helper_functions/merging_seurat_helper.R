library(tidyverse)
library(Seurat)
#library(DoubletFinder)

#' Read in all the filtered seurat objects to be merged together
#'
#' @param list of paths to the seurat objects
#'

read_seurat_objs <- function(paths_list){
	objs <- list()
	for(i in names(paths_list)){
		obj <- readRDS(paths_list[[i]])
		objs[[i]] <- RenameCells(obj, add.cell.id = i)
		#objs[[i]] <- obj
	}
	return(objs)
}

#' Summarize the filtering results by sample
#'
#' @param Filtered merged seurat object
#' @param Unfiltered merged seurat object
#'

summarize_filtering <- function(filtered_obj, unfiltered_obj){
	unfiltered_samp_info <- as.data.frame(table(unfiltered_obj@meta.data$orig.ident))
	colnames(unfiltered_samp_info) <- c("Sample", "Pre-Filtering")
	filtered_samp_info <- as.data.frame(table(filtered_obj@meta.data$orig.ident))
	colnames(filtered_samp_info) <- c("Sample", "Post-Filtering")
#	unfiltered_samp_info <- merge(unfiltered_samp_info, filtered_samp_info, by = "Sample")
	return(unfiltered_samp_info)
}


#' Generate violin plot template
#'
#' @param Merged seurat object
#' @param Feature to be plotted
#' @param Threshold value
#' @param y-axis title
#' 

generate_violin_plot <- function(merged_obj, feature, h, y_axis_title){
	h <- as.integer(h)
        plot1 <- VlnPlot(merged_obj, features = feature, pt.size=0, group.by="orig.ident", col=meta_data$sample_colour) +
                ggplot2::scale_y_continuous(name = y_axis_title, labels = scales::comma) +
                xlab(paste0(args$gr, " samples")) +
                ggtitle("") +
                theme_bw() +
                theme(legend.position = "none") +
                geom_hline(yintercept = h, linetype = 2, color = "red")
        return(plot1)
}


