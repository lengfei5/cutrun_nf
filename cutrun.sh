#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --qos=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-39
#SBATCH --output=/users/sean.montgomery/logs/cutrun_Tak1v5_%a.txt
# === end SBATCH directives ===

#Script to trim adapters then align CUT&RUN paired end data with Bowtie2
#Script adapted from Michael Borg

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
# sample_list=/users/sean.montgomery/inputfiles/cutrun_samplenames2-Tak1v3.txt
sample_list=/groups/berger/user/sean.montgomery/Documents/Scripts/input_files/cutrun_samplenames.yy.txt
# sample_list=$list_dir/cutrun_bamfiles.txt
#4. set directory containing sample folders
# work_folder=/scratch-cbe/users/sean.montgomery/cutrun/TakNv5
work_folder=/scratch-cbe/users/sean.montgomery/cutrun/Tak1v5
#5. bwt2 index for genome of interest
# index1=/groups/berger/lab/cluster_files/sean/TakNv5/TakNv5
index1=/groups/berger/lab/cluster_files/sean/Tak1v5/Tak1v5
#5. bwt2 index for spike-in genome
index2=/groups/berger/lab/cluster_files/NGS_annotations/backup_berger_common/hg19/hg19
#6. Data type "SE" or "PE"
datatype="PE"

F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
#4. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`

# Load the required modules
module load samtools/1.9-foss-2018b
module load bedtools/2.27.1-foss-2018b
module load cutadapt/1.18-foss-2018b-python-3.6.6
module load bowtie2/2.3.4.2-foss-2018b
module load deeptools/2.5.4-foss-2018b-python-2.7.15
module load picard/2.18.27-java-1.8
module load r/3.5.1-foss-2018b
# ... and then change to  working directory:
mkdir -p $work_folder/${NAME}
cd $work_folder
# === end ENVIRONMENT SETUP ===


# if [[ -f ${NAME}/${F} ]]; then
# 	echo "Bam file in folder"
# fi

## BAM to fastq
if [[ ${F##*.} == "bam" ]]; then
	cp /groups/berger/lab/Raw/demultiplexed/${F} $work_folder/${NAME}
        cp /groups/berger/lab/Raw/Ido_demult/HY7HLBGXC_20191204B_demux_1_M9395/${F} $work_folder/${NAME}
	# samtools sort -n -o ${NAME}.qsort -O bam ${F}
	# bedtools bamtofastq -i ${NAME}.qsort -fq ${NAME}.end1.fq -fq2 ${NAME}.end2.fq
	# echo "BAM converted to SAM for sample ... "
	# echo ${NAME}
	# rm ${NAME}/${NAME}.qsort
else
	gunzip ${F}_1.fastq.gz
	gunzip ${F}_2.fastq.gz
	mv ${F}_1.fastq ${NAME}.end1.fq
	mv ${F}_2.fastq ${NAME}.end2.fq
	echo "I tried my best with this SRA file"
fi

# if [[ -f ${F} ]]; then
# 	mv ${F} ${NAME}/
# fi


# #Pre-process files

# #NEBNext adaptor sequences
a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# # BAM file must be sorted/grouped by read name to keep records in the two output FASTQ files in the same order.
# # To sort the BAM file by query name use ---> samtools sort -n -o aln.qsort -O bam aln.bam

echo ".............................................................. "
echo "----> Converting BAM to fatsq files for sample ... " ${NAME}
samtools sort -n -o ${NAME}/${NAME}.qsort -O bam ${NAME}/${F}
if [ $datatype == "PE" ]; then
        bedtools bamtofastq -i ${NAME}/${NAME}.qsort -fq ${NAME}/${NAME}.end1.fq -fq2 ${NAME}/${NAME}.end2.fq
fi
if [ $datatype == "SE" ]; then
        bedtools bamtofastq -i ${NAME}/${NAME}.qsort -fq ${NAME}/${NAME}.fq
fi

# # cutadapt -a AACCGGTT -o output.fastq input.fastq --- removing sequences less than 1bp (self-ligated adapters)
echo ".............................................................. "
echo "---> Cutting adapters for sample ... " ${NAME}
if [ $datatype == "PE" ]; then
        cutadapt -a $a1 -A $a2 --minimum-length 1 -o ${NAME}/${NAME}.trim1.fq -p ${NAME}/${NAME}.trim2.fq ${NAME}/${NAME}.end1.fq ${NAME}/${NAME}.end2.fq > ${NAME}/${NAME}_cutadapt.txt
fi
if [ $datatype == "SE" ]; then
        cutadapt -a $a1 --minimum-length 1 -o ${NAME}/${NAME}.trim.fq  ${NAME}/${NAME}.fq > ${NAME}/${NAME}_cutadapt.txt
fi


# #End pre-processing files

# # align with bowtie2 fragments between 100bp and 700bp with 0 mismatches (-N)
echo ".............................................................. "
echo "---> Aligning sample ... " ${NAME}
if [ $datatype == "PE" ]; then
		bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700 -x $index1 -1 ${NAME}/${NAME}.trim1.fq -2 ${NAME}/${NAME}.trim2.fq -S ${NAME}/${NAME}.sam
		bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 700 -x $index2 -1 ${NAME}/${NAME}.trim1.fq -2 ${NAME}/${NAME}.trim2.fq -S ${NAME}/${NAME}.spike.sam
fi
if [ $datatype == "SE" ]; then
        bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700 -x $index1 -U ${NAME}/${NAME}.trim.fq -S ${NAME}/${NAME}.sam
	bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant -q --phred33 -I 10 -X 700 -x $index2 -U ${NAME}/${NAME}.trim.fq -S ${NAME}/${NAME}.spike.sam
fi

# # convert sam to bam
echo ".............................................................. "
echo "---> Converting SAM to BAM for sample ... " ${NAME}
samtools view -bS -o ${NAME}/${NAME}.aligned.bam ${NAME}/${NAME}.sam
samtools view -bS -o ${NAME}/${NAME}.spike.aligned.bam ${NAME}/${NAME}.spike.sam

# #Sorting and indexing aligned BAM file
echo ".............................................................. "
echo "---> Sorting by co-ordinates and indexing sample ... " ${NAME}.aligned.bam 
echo ${NAME}
samtools sort -o ${NAME}/${NAME}.aligned_sorted.bam ${NAME}/${NAME}.aligned.bam
samtools index ${NAME}/${NAME}.aligned_sorted.bam

samtools sort -o ${NAME}/${NAME}.spike.aligned_sorted.bam ${NAME}/${NAME}.spike.aligned.bam
samtools index ${NAME}/${NAME}.spike.aligned_sorted.bam

# #use samtools to filter away reads with MAPQ < 10
echo ".............................................................. "
echo "---> Filtering for mapQ > 10m sorting and indexing sample ... " ${NAME}.aligned_sorted.bam 
echo ${NAME}
samtools view -bq 10 ${NAME}/${NAME}.aligned_sorted.bam > ${NAME}/${NAME}.aligned_filtered.bam 
samtools sort -o ${NAME}/${NAME}.sorted.bam ${NAME}/${NAME}.aligned_filtered.bam 
samtools index ${NAME}/${NAME}.sorted.bam

samtools view -bq 10 ${NAME}/${NAME}.spike.aligned_sorted.bam > ${NAME}/${NAME}.spike.aligned_filtered.bam 
samtools sort -o ${NAME}/${NAME}.spike.sorted.bam ${NAME}/${NAME}.spike.aligned_filtered.bam 
samtools index ${NAME}/${NAME}.spike.sorted.bam

# #set tmp file for picard tools
export ptTMPDIR=$work_folder/picard_tmp
# # make sure the directory exists
mkdir -p $ptTMPDIR

# #removing duplicates and indexing resulting BAM file
echo ".............................................................. "
echo "---> Removing duplicates and indexing to generate ... " ${NAME}.sorted_uniq.bam
java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${NAME}/${NAME}.sorted.bam O=${NAME}/${NAME}.sorted_uniq.bam M=${NAME}/${NAME}.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
samtools index ${NAME}/${NAME}.sorted_uniq.bam

java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${NAME}/${NAME}.spike.sorted.bam O=${NAME}/${NAME}.spike.sorted_uniq.bam M=${NAME}/${NAME}.spike.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
samtools index ${NAME}/${NAME}.spike.sorted_uniq.bam

# #set tmp file for deeptools
export TMPDIR=$work_folder/deeptool_tmp
# # make sure the directory exists
mkdir -p $TMPDIR

# #normalise and convert BAM file to BIGWIG to represent coverage
echo ".............................................................. "
echo "---> Generating BW files for final BAM files with and without duplicates ... "
bamCoverage -b ${NAME}/${NAME}.sorted.bam -o ${NAME}/${NAME}.coverage_dups.bw --normalizeTo1x 239800000 --binSize=10
bamCoverage -b ${NAME}/${NAME}.sorted_uniq.bam -o ${NAME}/${NAME}.coverage_uniq.bw --normalizeTo1x 239800000 --binSize=10

bamCoverage -b ${NAME}/${NAME}.spike.sorted.bam -o ${NAME}/${NAME}.spike.coverage_dups.bw --normalizeTo1x 239800000 --binSize=10
bamCoverage -b ${NAME}/${NAME}.spike.sorted_uniq.bam -o ${NAME}/${NAME}.spike.coverage_uniq.bw --normalizeTo1x 239800000 --binSize=10

#convert BAM file to bedpe for further analysis
#set tmp file for sorting
export TMPDIR=$work_folder/tmp
# make sure the directory exists
mkdir -p $TMPDIR
if [ $datatype == "PE" ]; then
        echo ".............................................................. "
        echo "---> Generating BEDPE files for downstream analyses ... "
        samtools sort -n ${NAME}/${NAME}.sorted_uniq.bam | bamToBed -bedpe -i stdin | sort -k1,1V -k2,2n -T $TMPDIR > ${NAME}/${NAME}.sorted_uniq.bedpe
fi

# # generate alignment summary tables
# ### for aligned.bam (i,e to get raw mapping information)
echo ".............................................................. "
echo "---> Summarising alignment statistics for sample ... " ${NAME}
samtools stats ${NAME}/${NAME}.aligned_sorted.bam > ${NAME}/${NAME}.aligned_stats.txt
samtools stats ${NAME}/${NAME}.sorted.bam > ${NAME}/${NAME}.sorted_stats.txt
for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
        grep "^SN" ${NAME}/${NAME}.aligned_stats.txt| grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}/${NAME}.aligned_summary.txt
        grep "^SN" ${NAME}/${NAME}.sorted_stats.txt | grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}/${NAME}.sorted_summary.txt
done
rm ${NAME}/${NAME}.aligned_stats.txt
rm ${NAME}/${NAME}.sorted_stats.txt

samtools stats ${NAME}/${NAME}.spike.aligned_sorted.bam > ${NAME}/${NAME}.spike.aligned_stats.txt
samtools stats ${NAME}/${NAME}.spike.sorted.bam > ${NAME}/${NAME}.spike.sorted_stats.txt
for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
        grep "^SN" ${NAME}/${NAME}.spike.aligned_stats.txt| grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}/${NAME}.spike.aligned_summary.txt
        grep "^SN" ${NAME}/${NAME}.spike.sorted_stats.txt | grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}/${NAME}.spike.sorted_summary.txt
done
rm ${NAME}/${NAME}.spike.aligned_stats.txt
rm ${NAME}/${NAME}.spike.sorted_stats.txt

#filter fragments <150bp
samtools view -h ${NAME}/${NAME}.sorted_uniq.bam | awk '$9 > 150 || $9 < -150 || $1 ~ /^@/' | samtools view -bS > ${NAME}/${NAME}.sizednuc150.bam
samtools index ${NAME}/${NAME}.sizednuc150.bam

export TMPDIR=$work_folder/deeptool_tmp
bamCoverage -b ${NAME}/${NAME}.sizednuc150.bam -o ${NAME}/${NAME}.sizednuc150.bw --normalizeTo1x 239800000 --binSize=10
# bamCoverage -b ${NAME}/${NAME}.sizednuc150.bam -o ${NAME}/${NAME}.sizednuc150.bw --normalizeTo1x 119146348 --binSize=10

#convert BAM file to bedpe for further analysis
#set tmp file for sorting
export TMPDIR=$work_folder/tmp
# make sure the directory exists
mkdir -p $TMPDIR
if [ $datatype == "PE" ]; then
        echo ".............................................................. "
        echo "---> Generating BEDPE files for downstream analyses ... "
        samtools sort -n ${NAME}/${NAME}.sizednuc150.bam | bamToBed -bedpe -i stdin | sort -k1,1 -k2,2n -T $TMPDIR > ${NAME}/${NAME}.sizednuc150.bedpe
fi

# #Scale output according to spike-in
spike_frags=`samtools view ${NAME}/${NAME}.spike.sorted_uniq.bam | wc -l`
scale_factor=`echo "scale=5;10000 / $spike_frags" | bc`
samtools view -b ${NAME}/${NAME}.sizednuc150.bam | genomeCoverageBed -bg -scale $scale_factor -ibam stdin -g $chrom_sizes > ${NAME}/${NAME}.sizednuc150.norm.bed
