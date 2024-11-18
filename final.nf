#!/usr/bin/env nextflow

/*
 * Define default parameters
 */
params.genome  = "${projectDir}/data/ref_genome/ecoli_rel606.fasta"
params.reads   = "${projectDir}/data/trimmed_fastq_small/*_{1,2}_trim.fastq.gz"
params.results = "results"


process run_fastqc {
    publishDir "${params.results}/fastqc_reports", mode: 'symlink'
    input:
    tuple val(sampleId), path(reads)

    output:
    path "fastqc_reports"

    script:
    """
    mkdir -p fastqc_reports
    fastqc -o fastqc_reports ${reads}
    """
}

/*
 * Process 1A: Create FASTA genome index using SAMtools
 */
process index_genome_samtools {
    input:
    path genome

    output:
    path "${genome}.fai"

    script:
    """
    samtools faidx ${genome}
    """
}

/*
 * Process 1B: Create genome sequence dictionary using Picard
 */
process create_sequence_dict {
    input:
    path genome

    output:
    path "${genome.baseName}.dict"

    script:
    """
    picard CreateSequenceDictionary R=${genome} O=${genome.baseName}.dict
    """
}

/*
 * Process 1C: Create STAR genome index
 */
process index_genome_star {
    input:
    path genome

    output:
    path "star_index_dir"

    script:
    """
    mkdir -p star_index_dir

    STAR --runMode genomeGenerate \
         --genomeDir star_index_dir \
         --genomeFastaFiles ${genome} \
         --runThreadN ${task.cpus}
    """
}

/*
 * Process 2: Align RNA-Seq reads using STAR
 */
process align_reads_star {
    input:
    path genome
    path starIndexDir
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId),
          path('Aligned.sortedByCoord.out.bam'),
          path('Aligned.sortedByCoord.out.bam.bai')

    script:
    """
    # Align reads using STAR
    STAR --genomeDir ${starIndexDir} \
         --readFilesIn ${reads} \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:${sampleId} LB:library PL:illumina PU:machine SM:sample

    # Index the sorted BAM file
    samtools index Aligned.sortedByCoord.out.bam
    """
}

/*
 * Process 3: Split reads on 'N' in CIGAR string using GATK
 */
process split_reads_gatk {
    input:
    path genome
    path genomeIndex
    path genomeDict
    tuple val(sampleId), path(bam), path(bai)

    output:
    tuple val(sampleId), path('split.bam'), path('split.bai')

    script:
    """
    java -jar /usr/gitc/GATK35.jar -T SplitNCigarReads \
                                   -R ${genome} -I ${bam} \
                                   -o split.bam \
                                   -rf ReassignOneMappingQuality \
                                   -RMQF 255 -RMQT 60 \
                                   -U ALLOW_N_CIGAR_READS
    """
}

/*
 * Process 4: Variant Calling using GATK HaplotypeCaller
 */
process call_variants_gatk {
    publishDir "${params.results}/variant_calls", mode: 'symlink'
    input:
    path genome
    path genomeIndex
    path genomeDict
    tuple val(sampleId), path(bam), path(bai)

    output:
    tuple val(sampleId), path('final.vcf')

    script:
    """
    echo "${bam.join('\n')}" > bam.list

    # Variant calling with HaplotypeCaller
    java -jar /usr/gitc/GATK35.jar -T HaplotypeCaller \
                                   -R ${genome} -I bam.list \
                                   -dontUseSoftClippedBases \
                                   -stand_call_conf 20.0 \
                                   -o output.gatk.vcf.gz

    # Variant filtering with GATK
    java -jar /usr/gitc/GATK35.jar -T VariantFiltration \
                                   -R ${genome} -V output.gatk.vcf.gz \
                                   -window 35 -cluster 3 \
                                   -filterName FS -filter "FS > 30.0" \
                                   -filterName QD -filter "QD < 2.0" \
                                   -o final.vcf
    """
}

/*
 * Process 5: Filter VCF for quality metrics
 */
process filter_vcf {
    publishDir "${params.results}/filtered_vcf", mode: 'symlink'
    input:
    tuple val(sampleId), path('final.vcf')

    output:
    tuple val(sampleId), path('filtered.vcf')

    script:
    """
    # Filter VCF based on depth (DP >= 8)
    grep -v '#' final.vcf | awk '\$7~/PASS/' | perl -ne 'chomp(\$_); \
        (\$dp)=\$_=~/DP\\=(\\d+)\\;/; if(\$dp>=8){print \$_."\\n"};' > filtered.vcf
    """
}

/*
 * Process 6: Coverage Analysis using SAMtools
 */
process coverage_analysis {
    publishDir "${params.results}/coverage_reports", mode: 'symlink'
    input:
    tuple val(sampleId), path(bam)

    output:
    path("${sampleId}_coverage.txt")

    script:
    """
    samtools depth ${bam} > ${sampleId}_coverage.txt
    """
}

/*
 * Workflow: Define the pipeline execution steps
 */
workflow {

    // Channel for paired-end reads
    reads_ch = Channel.fromFilePairs(params.reads)
    run_fastqc(reads_ch)
    // Step 1: Index the genome
    index_genome_samtools(params.genome)
    create_sequence_dict(params.genome)
    index_genome_star(params.genome)

    // Step 2: Align reads to the genome
    align_reads_star(params.genome, index_genome_star.out, reads_ch)

    // Step 3: Split reads on 'N' in CIGAR string
    split_reads_gatk(params.genome,
                     index_genome_samtools.out,
                     create_sequence_dict.out,
                     align_reads_star.out)

    // Group BAM files by sample ID
    split_reads_gatk.out
        .map { sampleId, bam, bai -> tuple(sampleId, bam, bai) }
        .groupTuple()
        .set { grouped_bams }

    // Step 4: Variant calling
    call_variants_gatk(params.genome,
                       index_genome_samtools.out,
                       create_sequence_dict.out,
                       grouped_bams)

    // Step 5: Filter VCF file for quality
    filter_vcf(call_variants_gatk.out)

    // Coverage analysis
    grouped_bams.map { sampleId, bam, _ -> tuple(sampleId, bam) }
                .set { bam_files }
    coverage_analysis(bam_files)
}