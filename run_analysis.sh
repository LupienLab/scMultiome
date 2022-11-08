#!/bin/bash

module load R/4.1.0

# Create some output directories
mkdir -p ./logs
mkdir -p ./seurat_objects/unfiltered_seurat_obj/merged_seurat_obj
mkdir -p ./seurat_objects/filtered_seurat_obj/merged_seurat_obj
mkdir -p ./plots/QC

# Arguments:
## $1: Path to the submit samples csv file
## $2: An array of the run steps - options: ("obj_gen" "merge")
## $3: Working with unfiltered or filtered seurat objects? Applicable if $2 = ("merge")

file=$1
step=$2
obj=$3

reads_cutoffs=(2000 1500 1000 500 250)

# Create unfiltered, run filtering QCs and generate filtered seurat objects
if [ $step = "obj_gen" ]
then
{
read
while IFS="," read -r fastq_name sample_ID group sample_colour mito_threshold HB_threshold nFeature_threshold nCount_threshold
do
   mkdir -p ./plots/QC/$sample_ID

   echo "Submitting job for" $sample_ID

   sbatch -p all --mem=30G -o ./logs/$sample_ID'.out' -e ./logs/$sample_ID'.err' -t 7:0:0 -J $sample_ID'_seurat' --wrap "Rscript 01-Generating_Seurat_Obj.R -m /cluster/home/sbahl/Catherine_DIR/metrib/data/$group'_data/'$fastq_name -s $sample_ID -col $sample_colour -mito $mito_threshold -HB $HB_threshold -gene $nFeature_threshold -count $nCount_threshold > ./logs/01-Generating_Seurat_Obj_$sample_ID'.Rout'"  

done 
} < $file
fi

# Merge the seurat objects according to the groups in the sample sheet
if [ $step = "merge" ]
then
	groups=$(awk -F ',' '(NR>1) {print $3}' $file | sort -u)

	for group in ${groups[@]}
	do
		mkdir -p ./plots/QC/$group
		echo "Running merging for "$group" samples"
		sample_array=$(awk -F ',' -v gr="$group" '$3 ==gr { print $2 }' $file)
		echo "Merging these samples" ${samples[@]}

		sbatch -p all --mem=30G -o ./logs/$group'_merge.out' -e ./logs/$group'_merge.err' -t 1:0:0 -J $group'_merge' --wrap "Rscript 02-Merging_Seurat_Obj.R -f $file -gr $group -obj $obj > ./logs/02-Merging_Seurat_Obj_$group'.Rout'" 

		#for cutoff in ${reads_cutoffs[@]}
		#do
	#		echo "Running with reads cutoff" $cutoff
	#		sbatch -p all --mem=30G -o ./logs/$group'_'$cutoff'_merge.out' -e ./logs/$group'_'$cutoff'_merge.err' -t 1:0:0 -J $group'_'$cutoff'_merge' --wrap "Rscript 03-QC_plots_merged_obj.R -f $file -r $cutoff -gr $group -obj $obj > ./logs/03-QC_plots_merged_obj_$group'_'$cutoff'.Rout'"
			#Rscript 03-QC_plots_merged_obj.R -f $file -gr $group -r $cutoff -obj $obj > ./logs/03-QC_plots_merged_obj_$group'_'$cutoff'.Rout'
	#	done
		#sbatch -p all --mem=30G -o ./logs/$group'_merge.out' -e ./logs/$group'_merge.err' -t 1:0:0 -J $group'_merge' --wrap "Rscript 03-QC_plots_merged_obj.R -f $file -gr $group -obj $obj > ./logs/03-QC_plots_merged_obj_$group'.Rout'"

	done
fi

quit()

##### Debugging below #####



