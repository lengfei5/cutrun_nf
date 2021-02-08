#!/usr/bin/env nextflow

/*
*****************************
 * Pipeline - smallRNA-meth *
 *                          *
 * Thomas R Burkard         *
 * IMP/IMBA Bioinformatics  *
 ****************************
*/


log.info """\
         =============================
         Pipeline - smallRNA-meth
         =============================

         outdir: ${params.outdir}
         reads: ${params.reads}
         genome: ${params.genome}
         adapter1: ${params.adapter1}
         adapter2: ${params.adapter2}
         tss: ${params.tss}

         """
         .stripIndent()

/*
 * Input parameters validation
 */

if (params.genome)
{
  genome_file = file(params.genome)
  //if( !genome_file.exists() ) exit 1, "Genome file doesn't exist: ${genome_file}"

} else {
  exit 1, "Missing genome file"
}

if (params.tss)
{
  tss_file = file(params.tss)
  if( !tss_file.exists() ) exit 1, "tss file doesn't exist: ${tss_file}"
} else {
  tss_file = file("NA")
}

/*
 * Validate input files
 */

/*
 * Create a channel for read files
 */
Channel
    .fromPath( params.reads)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { file -> tuple(file.baseName, file) }
    .set { read_files }


/*
 * bam to fastq and pair splitting
 */
process splitfastq {
    tag "Channel: ${name}"

    publishDir "${params.outdir}/FASTQs", mode: 'copy', pattern: '*.fastq'

    input:
        set val(name), file(bam) from read_files

    output:
        set name, file("*.fastq") into splited_fastq
        set name, file("cntTotal.txt") into cnt_total

    script:
    """
        ml load bedtools/2.25.0-foss-2018b
        ml load samtools/1.10-foss-2018b

        samtools view -c ${bam} > cntTotal.txt

        bamToFastq -i ${bam} -fq ${name}.fastq
        bash ${baseDir}/scripts/deinterleave_fastq.sh < ${name}.fastq ${name}_R1.fastq ${name}_R2.fastq

        rm ${name}.fastq

    """
}


/*
 * cut adapter
 */
process  cutadapt {
    tag "Channel: ${name}"

    publishDir "${params.outdir}/cutadapt", mode: 'copy', pattern: '*.err'

    input:
        set val(name), file(fastq) from splited_fastq

    output:
        set name, file("${name}_R*.trim.fastq") into fastq_trimmed
        set name, file("cutadapt.${name}.err") into stat_cutadapt

    script:
    """
        module load python/2.7.15-gcccore-7.3.0-bare
        module load cutadapt/1.18-foss-2018b-python-2.7.15

        cutadapt --minimum-length ${params.minLength} --overlap ${params.overlapLength} -a ${params.adapter1} -A ${params.adapter2} \
        -o ${name}_R1.trim.fastq -p ${name}_R2.trim.fastq ${name}_R1.fastq ${name}_R2.fastq > cutadapt.${name}.err

    """
}

/*
 * Align with bowite2 with some particular parameters for cut_and_run
 */
process align {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/aligned_bams", mode: 'copy'

    input:
    set name, file(fastq) from fastq_trimmed

    output:
      set name, file("*.bam") into aligned_bam, aligned_bam2
      file "${name}.bam*"

    script:
    """
    ml load bowtie2/2.3.5.1-foss-2018b
    ml load samtools/1.10-foss-2018b
    bowtie2 --end-to-end --very-sensitive --phred33 -I 10 -q --no-mixed -X 700 --dovetail --no-discordant -x ${params.genome} -1 ${name}_R1.trim.fastq -2 ${name}_R2.trim.fastq | \
    samtools view -bSu - | samtools sort - -o ${name}.bam

    samtools index ${name}.bam

    """
}


/*
 * filter algined bam, duplicated removal
 */
process filter_rmdup_bam {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/filtered_bams", mode: 'copy'

    input:
    set name, file(bam) from aligned_bam

    output:
      set name, file("${name}_uniq.rmdup.bam") into filterd_bam

      set name, file ("${name}_cnt_trimmed.txt") into trimmed_stat
      set name, file ("${name}_cnt_mapped.txt") into mapped_stat
      set name, file ("${name}_cnt_uniq.txt") into unique_stat
      set name, file ("${name}_cnt_uniq.rmdup.txt") into usable_stat

      file "${name}_uniq.rmdup.bam"
      file '*.pdf'

    script:
    """
      ml load samtools/1.10-foss-2018b
      ml load picard/2.20.6--0-biocontainers

      samtools view -h -q 30 ${bam} > ${name}.unsorted.bam
      samtools sort -o ${name}_uniq.bam ${name}.unsorted.bam
      samtools index ${name}_uniq.bam
      rm ${name}.unsorted.bam

      picard MarkDuplicates INPUT=${name}_uniq.bam OUTPUT=${name}_uniq.rmdup.unsorted.bam METRICS_FILE=${name}_picard.rmDup.txt \
      ASSUME_SORTED=true REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.2

      samtools sort -o ${name}_uniq.rmdup.bam ${name}_uniq.rmdup.unsorted.bam
      samtools index ${name}_uniq.rmdup.bam
      rm ${name}_uniq.rmdup.unsorted.bam

      picard CollectInsertSizeMetrics HISTOGRAM_FILE=${name}_fragsize.pdf OUTPUT=${name}_frag_size.txt \\
      METRIC_ACCUMULATION_LEVEL=ALL_READS INCLUDE_DUPLICATES=false INPUT=${name}_uniq.rmdup.bam

      samtools view -c ${bam} > ${name}_cnt_trimmed.txt
      samtools view -c -F 4 ${bam} > ${name}_cnt_mapped.txt
      samtools view -c ${name}_uniq.bam > ${name}_cnt_uniq.txt
      samtools view -c ${name}_uniq.rmdup.bam > ${name}_cnt_uniq.rmdup.txt

    """
}


/*
 * make bigwig with deeptools
 */
process bigwig_deeptools {

    tag "Channel: ${name}"
    publishDir "${params.outdir}/bigwigs", mode: 'copy'

    input:
    set name, file(bam) from aligned_bam2

    output:
      file "*.bw"

    shell:
    """
    ml load samtools/1.10-foss-2018b
    ml load deeptools/3.1.2-foss-2018b-python-2.7.15

    samtools index -c -m 14 !{bam}

    samtools view -h !{bam} | awk '\$9 > 150 || \$9 < -150 || \$1 ~ /^@/' | samtools view -bS > !{name}.sizednuc150.bam
    samtools index !{name}.sizednuc150.bam

    bamCoverage -b !{bam} -o !{name}.bw --outFileFormat=bigwig --normalizeUsing CPM --ignoreDuplicates --binSize 20
    bamCoverage -b !{name}.sizednuc150.bam -o !{name}.sizednuc150.bw --outFileFormat=bigwig --normalizeUsing CPM --ignoreDuplicates --binSize 20

    """
}

/*
 * make statTable
 */
 cnt_total.concat(trimmed_stat, mapped_stat, unique_stat, usable_stat)
     .groupTuple()
     .map{ stat1, stat2 -> [stat1, stat2[0], stat2[1], stat2[2], stat2[3], stat2[4]] }
     .set{ cntStat_files }

process statTable {

    tag "Channel: ${name}"

    publishDir "${params.outdir}/countStat", mode: 'copy'
    input:
        set name, file(cnt_total), file(trimmed), file(mapped), file(unique_stat), file(usable_stat) from cntStat_files

    output:
        file "${name}.countStat.txt" into cnt_stat

    script:
    """
    echo -e "Name\tTotal\ttrimmed\taligned\tunique\tunique.rmdup" > ${name}.countStat.txt
    TOTAL=`cat ${cnt_total}`
    TRIMMED=`cat ${name}_cnt_trimmed.txt`
    MAPPED=`cat ${name}_cnt_mapped.txt`
    UNIQUE=`cat ${name}_cnt_uniq.txt`
    USABLE=`cat ${name}_cnt_uniq.rmdup.txt`

    echo -e "${name}\t\$TOTAL\t\$TRIMMED\t\$MAPPED\t\$UNIQUE\t\$USABLE" >> ${name}.countStat.txt

    """
}

/*
 *  stat table of all samples
 */
process countTable {

    publishDir "${params.outdir}/result", mode: 'copy'

    input:
	      file "countStat/*" from cnt_stat.collect()

    output:
	      file "countStatTable.txt"

    script:
    """
    ml load r/3.6.2-foss-2018b
    cp $baseDir/scripts/countTable.R ./countTable.R
    Rscript countTable.R

    """
}



workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Execution status      : ${ workflow.success ? 'succeeded' : 'failed' }"
    println "Duration : ${workflow.duration}"
}
