#!/usr/bin/env nextflow

/*
 * Define the default parameters
 */
params.genome     = "${projectDir}/data/ref_genome/ecoli_rel606.fasta"
params.reads      = "${projectDir}/data/trimmed_fastq_small/*_{1,2}_trim.fastq.gz"
params.results    = "results"

/*
 * Process 1A: Create a FASTA genome index with samtools
 */
process prepare_genome_samtools {
    container 'quay.io/biocontainers/samtools:1.3.1--h0cf4675_11'

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
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process prepare_genome_picard {
    container 'quay.io/biocontainers/picard:1.141--hdfd78af_6'

    input:
    path genome

    output:
    path "${genome.baseName}.dict"

    script:
    """
    picard CreateSequenceDictionary R= ${genome} O= ${genome.baseName}.dict
    """
}

/*
* Process 1C: Create the genome index file for STAR
*/
process prepare_star_genome_index {
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

    input:
    path genome

    output:
    path 'genome_dir'

    script:
    """
    mkdir -p genome_dir

    STAR --runMode genomeGenerate \
         --genomeDir genome_dir \
         --genomeFastaFiles ${genome} \
         --runThreadN ${task.cpus}
    """
}

/*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process rnaseq_mapping_star {
    container 'quay.io/biocontainers/mulled-v2-52f8f283e3c401243cee4ee45f80122fbf6df3bb:e3bc54570927dc255f0e580cba1789b64690d611-0'

    input:
    path genome
    path genomeDir
    tuple val(replicateId), path(reads)

    output:
    tuple val(replicateId),
          path('Aligned.sortedByCoord.out.bam'),
          path('Aligned.sortedByCoord.out.bam.bai')

    script:
    """
    # ngs-nf-dev Align reads to genome
    STAR --genomeDir ${genomeDir} \
         --readFilesIn ${reads} \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999

    # 2nd pass (improve alignments using table of splice
    # junctions and create a new index)
    mkdir -p genomeDir
    STAR --runMode genomeGenerate \
         --genomeDir genomeDir \
         --genomeFastaFiles ${genome} \
         --sjdbFileChrStartEnd SJ.out.tab \
         --sjdbOverhang 75 \
         --runThreadN ${task.cpus}

    # Final read alignments
    STAR --genomeDir genomeDir \
         --readFilesIn ${reads} \
         --runThreadN ${task.cpus} \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:${replicateId} LB:library PL:illumina \
                            PU:machine SM:GM12878

    # Index the BAM file
    samtools index Aligned.sortedByCoord.out.bam
    """
}

/*
 * Process 3: GATK Split on N
 */

process rnaseq_gatk_splitNcigar {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${replicateId}"

    input:
    path genome
    path index
    path genome_dict
    tuple val(replicateId),
          path(bam),
          path(bai)

    output:
    tuple val(replicateId), path('split.bam'), path('split.bai')

    script:
    """
    # SplitNCigarReads and reassign mapping qualities
    java -jar /usr/gitc/GATK35.jar -T SplitNCigarReads \
                                   -R ${genome} -I ${bam} \
                                   -o split.bam \
                                   -rf ReassignOneMappingQuality \
                                   -RMQF 255 -RMQT 60 \
                                   -U ALLOW_N_CIGAR_READS
    """
}

/*
 * Process 5: GATK Variant Calling
 */

process rnaseq_call_variants {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${sampleId}"

    input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(bam), path(bai)

    output:
    tuple val(sampleId), path('final.vcf')

    script:
    """
    echo "${bam.join('\n')}" > bam.list

    # Variant calling
    java -jar /usr/gitc/GATK35.jar -T HaplotypeCaller \
                    -R ${genome} -I bam.list \
                    -dontUseSoftClippedBases \
                    -stand_call_conf 20.0 \
                    -o output.gatk.vcf.gz

    # Variant filtering
    java -jar /usr/gitc/GATK35.jar -T VariantFiltration \
                                   -R ${genome} -V output.gatk.vcf.gz \
                                   -window 35 -cluster 3 \
                                   -filterName FS -filter "FS > 30.0" \
                                   -filterName QD -filter "QD < 2.0" \
                                   -o final.vcf
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)

    prepare_genome_samtools(params.genome)
    prepare_genome_picard(params.genome)
    prepare_star_genome_index(params.genome)

    rnaseq_mapping_star(params.genome, prepare_star_genome_index.out, reads_ch)

    rnaseq_gatk_splitNcigar(params.genome,
                            prepare_genome_samtools.out,
                            prepare_genome_picard.out,
                            rnaseq_mapping_star.out)

    // Group BAMs by sampleId after extracting it from replicateId
    rnaseq_gatk_splitNcigar.out
        .map { replicateId, bam, bai ->
            sampleId = replicateId.replaceAll(/[12]$/,'')
            tuple(sampleId, bam, bai)
        }
        .groupTuple()
        .set { split_bams }

    rnaseq_call_variants(params.genome,
                         prepare_genome_samtools.out,
                         prepare_genome_picard.out,
                         split_bams)
}