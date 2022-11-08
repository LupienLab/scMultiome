# This script is to merge all the seurat objects in respective groups

library(tidyverse)
library(Seurat)
library(argparse)

source("./helper_functions/merging_seurat_helper.R")

#### Collect all the arguments ####################################################################

parser <- ArgumentParser()

parser$add_argument("-f", "--meta_data", required=TRUE, help="The path to the .csv meta-data file")
parser$add_argument("-gr", "--group", required=TRUE, help="Group name for samples to be merged")
parser$add_argument("-obj", "--obj_processed", required=TRUE, help="Should the 'filtered' or 'unfiltered' object be merged")

print(parser)

args <- parser$parse_args()
print(args)

###################################################################################################

meta <- read.table(args$meta_data, header = TRUE, sep = ',')
meta <- meta %>% filter(group == args$group)
print(head(meta))

# Read in all the filtered seurat objects
sample_objs <- list()

for(sample in meta$sample_ID){
	sample_obj_path <- paste0("./seurat_objects/", args$obj_processed, "_seurat_obj/", sample, "_", args$obj_processed, "_seurat_obj.rds")
	sample_objs[[sample]] <- sample_obj_path
}

print(sample_objs)

RNA_objs <- read_seurat_objs(sample_objs)
print(RNA_objs)

print(head(rownames(RNA_objs[[names(RNA_objs)[1]]]@meta.data)))

merge_obj <- reduce(RNA_objs, merge)
print('Done merging combined seurat matrix')
print(merge_obj)

saveRDS(merge_obj, paste0("./seurat_objects/", args$obj_processed, "_seurat_obj/merged_seurat_obj/", args$gr, "_", args$obj_processed, "_merged_seurat_obj.rds"))

if(args$obj_processed == "filtered"){
        print("Can summarize the number of cells pre-doublets")
	unfiltered_obj <- readRDS(paste0("./seurat_objects/unfiltered_seurat_obj/", sample, "_unfiltered_seurat_obj.rds"))
	print(table(unfiltered_obj@meta.data$orig.ident))
	print(table(merge_obj@meta.data$orig.ident))

        total_cells <- summarize_filtering(merge_obj, unfiltered_obj)
	print(total_cells)

	write.table(total_cells, paste0("./plots/QC/", args$gr, "/", args$gr, "_filtering_summary.tsv"), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
} else {print("Now add in filtering parameters to the csv meta data file")}


sessionInfo()

