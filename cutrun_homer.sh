#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-39
#SBATCH --output=/users/sean.montgomery/logs/cutrun_homer_Tak1v5_%a.txt
# === end SBATCH directives ===

#Script to trim adapters then align CUT&RUN paired end data with Bowtie2
#Script adapted from Michael Borg

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/cutrun_samplenames.yy.txt
#4. set directory containing sample folders
work_folder=/scratch-cbe/users/sean.montgomery/cutrun/Tak1v5
#6. Data type "SE" or "PE"
datatype="PE"

F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
#4. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`

# Load the required modules
module load samtools/1.9-foss-2018b
module load homer/4.10-foss-2018b
# ... and then change to  working directory:
cd $work_folder
# === end ENVIRONMENT SETUP ===

mkdir -p ${NAME}/homer
export TMPDIR=$work_folder/tmp
mkdir -p $TMPDIR

if [[ -f ${F} ]]; then
	mv ${F} ${NAME}/
fi

# Make tag directory
if [ $datatype == "PE" ]; then
        # makeTagDirectory ${NAME}/homer/ -sspe -single ${NAME}/${NAME}.sorted_uniq.bam
        makeTagDirectory ${NAME}/homer/ -sspe -single ${NAME}/${NAME}.sizednuc150.bam
fi
if [ $datatype == "SE" ]; then
        makeTagDirectory ${NAME}/homer/ -single ${NAME}/${NAME}.sizednuc150.bam
fi

# Visualizing experiments with a genome browser

# Find enriched peaks and regions
# Play with these parameters to find something that fits best; may need to vary per mark
# findPeaks ${NAME}/homer/ -style histone -size 250 -minDist 500 -o ${NAME}/homer/${NAME}.regions.narrow.txt
findPeaks ${NAME}/homer/ -style histone -size 750 -minDist 1500 -o ${NAME}/homer/${NAME}.regions.broad.txt


# # Convert output regions to viewable format
# pos2bed.pl ${NAME}/homer/${NAME}.regions.narrow.txt > ${NAME}/homer/${NAME}.regions.narrow.bed
# grep -v '#' ${NAME}/homer/${NAME}.regions.narrow.bed | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/${NAME}.regions.narrow.bedgraph
pos2bed.pl ${NAME}/homer/${NAME}.regions.broad.txt > ${NAME}/homer/${NAME}.regions.broad.bed
grep -v '#' ${NAME}/homer/${NAME}.regions.broad.bed | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/${NAME}.regions.broad.bedgraph

# Run locally after moving files
# for f in *bedgraph; do NAME=`basename $f .bedgraph`; bedGraphToBigWig $f ~/Documents/HiC/HiC-Pro_2.9.0/annotation/MpTak1v3.chrom.sizes bw/$NAME.bw; done
# for f in *bedgraph; do NAME=`basename $f .bedgraph`; bedGraphToBigWig $f ~/Documents/HiC/HiC-Pro_2.9.0/annotation/Tak1v4.chrom.sizes bw/$NAME.bw; done


# # Make tag directory
# if [ $datatype == "PE" ]; then
#         makeTagDirectory ${NAME}/homer/all -sspe -single ${NAME}/${NAME}.sorted_uniq.bam
#         makeTagDirectory ${NAME}/homer/tf -sspe -single ${NAME}/${NAME}.tf150.bam
# fi


# # Visualizing experiments with a genome browser

# # Find enriched peaks and regions
# # Play with these parameters to find something that fits best; may need to vary per mark
# findPeaks ${NAME}/homer/all -style histone -size 250 -minDist 500 -o ${NAME}/homer/all/${NAME}.regions.narrow.all.txt
# findPeaks ${NAME}/homer/all -style histone -size 750 -minDist 1500 -o ${NAME}/homer/all/${NAME}.regions.broad.all.txt
# findPeaks ${NAME}/homer/tf -style histone -size 250 -minDist 500 -o ${NAME}/homer/tf/${NAME}.regions.narrow.tf.txt
# findPeaks ${NAME}/homer/tf -style histone -size 750 -minDist 1500 -o ${NAME}/homer/tf/${NAME}.regions.broad.tf.txt

# # # Convert output regions to viewable format
# pos2bed.pl ${NAME}/homer/all/${NAME}.regions.narrow.all.txt > ${NAME}/homer/all/${NAME}.regions.narrow.all.bed
# grep -v '#' ${NAME}/homer/all/${NAME}.regions.narrow.all.bed | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/all/${NAME}.regions.narrow.all.bedgraph
# pos2bed.pl ${NAME}/homer/all/${NAME}.regions.broad.all.txt > ${NAME}/homer/all/${NAME}.regions.broad.all.bed
# grep -v '#' ${NAME}/homer/all/${NAME}.regions.broad.all.bed | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/all/${NAME}.regions.broad.all.bedgraph

# pos2bed.pl ${NAME}/homer/tf/${NAME}.regions.narrow.tf.txt > ${NAME}/homer/tf/${NAME}.regions.narrow.tf.bed
# grep -v '#' ${NAME}/homer/tf/${NAME}.regions.narrow.tf.bed | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/tf/${NAME}.regions.narrow.tf.bedgraph
# pos2bed.pl ${NAME}/homer/tf/${NAME}.regions.broad.tf.txt > ${NAME}/homer/tf/${NAME}.regions.broad.tf.bed
# grep -v '#' ${NAME}/homer/tf/${NAME}.regions.broad.tf.bed | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/tf/${NAME}.regions.broad.tf.bedgraph

# ##Create bedgraph files with peak scores
# grep -v '#' ${NAME}/homer/${NAME}.regions.narrow.bed > ${NAME}/homer/${NAME}.regions.narrow.2.bed
# grep -v '#' ${NAME}/homer/${NAME}.regions.narrow.txt | awk '{print $8}' > ${NAME}/homer/${NAME}.regions.narrow.scores.txt
# paste ${NAME}/homer/${NAME}.regions.narrow.2.bed ${NAME}/homer/${NAME}.regions.narrow.scores.txt | awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$7}' | LC_COLLATE=C sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/homer/${NAME}.regions.narrow.scores.bedgraph
# rm ${NAME}/homer/${NAME}.regions.narrow.2.bed
# rm ${NAME}/homer/${NAME}.regions.narrow.scores.txt

