# Variant Calling Analysis Pipeline

## Introduction

This RNA-Seq analysis pipeline automates the key stages of processing high-throughput sequencing data, from quality control of raw reads to variant calling and coverage analysis. RNA-Seq is widely used to investigate gene expression and detect genetic variants, but the analysis can be complex and involve several interdependent tools. This pipeline integrates commonly used bioinformatics software tools, including FastQC, SAMtools, Picard, STAR, and GATK, to streamline the entire workflow.

By leveraging Docker containers, the pipeline provides a consistent environment for each tool, ensuring reproducible results and simplifying the management of software dependencies. The primary aim of this pipeline is to provide an efficient, reproducible method for processing RNA-Seq data, handling tasks such as read alignment, genome indexing, and variant calling. The integration of Docker ensures that each tool operates in a stable environment, reducing compatibility issues across different systems. This pipeline is particularly suitable for projects requiring consistent processing of multiple RNA-Seq samples, enabling researchers to efficiently obtain high-quality outputs for downstream analyses.

This script was developed and tested in a Gitpod environment, an online IDE that offers a pre-configured and containerized workspace. Using Gitpod allows for seamless development with all dependencies installed and configured. Users can easily open the Gitpod workspace directly from the GitHub repository by clicking the button on the repository page.

## Features

- **Automated Workflow**: Streamlines RNA-Seq data processing from raw reads to variant calling.
- **Docker Integration**: Ensures consistent environments for each tool, enhancing reproducibility.
- **Multi-tool Support**: Integrates FastQC, SAMtools, Picard, STAR, and GATK.
- **Scalability**: Suitable for processing multiple RNA-Seq samples consistently.
- **Comprehensive Outputs**: Generates quality control reports, variant calls, filtered VCFs, and coverage reports.

---

## Usage

### Prerequisites

Installing this pipeline requires the following software:

- **Git**: For cloning the repository. [Installation Guide](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- **Nextflow**: Workflow management system. [Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- **Docker**: For containerized execution. [Installation Guide](https://docs.docker.com/get-docker/)

### Installation

1. **Clone the Repository**  
   Open a terminal and run the following command:

   ```bash
   git clone https://github.com/yiyang2002/501.git](https://github.com/yiyang2002/variant_calling.git
   ```
   
2. **Navigate to the Pipeline Directory**
   
   ```bash
   cd variant_calling
   ```

## Docker Images

The pipeline uses Docker containers for each process, as specified in nextflow.config. This ensures a consistent environment and eliminates issues related to software dependencies and version conflicts. The following Docker images are used:
	
•	FastQC: biocontainers/fastqc:v0.11.9_cv8 [1]
 
•	SAMtools: quay.io/biocontainers/samtools:1.3.1--h0cf4675_11 [2]
 
•	Picard: quay.io/biocontainers/picard:1.141--hdfd78af_6 [3]
 
•	STAR: quay.io/biocontainers/star:2.7.10b--h6b7c446_1 [4]
 
•	GATK: quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571 [5]
 

## Running the Pipeline

Execute the pipeline using Nextflow:

nextflow run final.nf

	•	Note: Ensure that Docker is running on your system before executing the command.
	•	The pipeline will take around 3 minutes to complete.

## Gitpod Environment

This pipeline was developed and tested in Gitpod, which provides a pre-configured development environment.

Launching Gitpod Workspace (Optional)

  1.	Open the GitHub Repository
Navigate to the project’s GitHub repository.

  2.	Launch Gitpod
Click on the Gitpod button on the repository page.

  3.	Run the Pipeline
Inside Gitpod, execute the pipeline as described in the Running the Pipeline section.

---

Pipeline Overview
The following is a simplified representation of the workflow in a Directed Acyclic Graph (DAG):

![image](https://github.com/user-attachments/assets/4e7af725-abb0-4885-80a6-4932bbafa47b)


---

## Input Data

The input data for this RNA-Seq pipeline consists of paired-end FASTQ files and a reference genome, originating from a long-term evolution experiment (LTEE) with Escherichia coli. Specifically, the dataset features a clone from the 50,000th generation, submitted by The University of Texas at Austin as part of a study on E. coli genome evolution over 50,000 generations.
### Note:
Both the paired-end FASTQ files and the reference genome (*E. coli REL606*) are already included in the repository under the `data/` directory. If needed, you can manually follow the steps below to regenerate the files.

### Subsampling Command

To reduce the dataset size for efficient analysis, a subset of 100,000 reads was extracted from each FASTQ file using seqtk:

```bash
$ seqtk sample -s100 SRR2584866_1.fastq.gz 100000 > SRR2584866_1_trim.fastq
$ seqtk sample -s100 SRR2584866_2.fastq.gz 100000 > SRR2584866_2_trim.fastq
```

### Download the reference genome for E. coli REL606 
```bash
$ mkdir -p data/ref_genome
$ curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
$ gunzip data/ref_genome/ecoli_rel606.fasta.gz
```
---

## Output Data

### Results Output

These are the final outputs generated by the pipeline, stored in the `results/` directory. Each subdirectory contains specific types of results critical for downstream analysis.

---

### 1. **Coverage Reports**
   - **Location**: `results/coverage_reports/`
   - **File**: `SRR2584866_coverage.txt`
   - **Description**: 
     - This text file provides detailed information about sequencing depth across the reference genome.
     - Includes per-base coverage values, which help assess whether there is sufficient coverage for reliable variant detection and expression analysis.

### 2. **Quality Control Reports**
   - **Location**: `results/fastqc_reports/`
   - **Files**:
     - `SRR2584866_1_trim_fastqc.html`: A comprehensive HTML report summarizing quality metrics for the first paired-end read file (`SRR2584866_1_trim.fastq.gz`).
     - `SRR2584866_1_trim_fastqc.zip`: Compressed file containing additional data, including raw statistics and visualizations for the first paired-end read.
     - `SRR2584866_2_trim_fastqc.html`: Similar to the above, but for the second paired-end read file (`SRR2584866_2_trim.fastq.gz`).
     - `SRR2584866_2_trim_fastqc.zip`: Compressed data for the second paired-end read.
   - **Description**:
     - These reports include critical metrics such as per-base sequence quality, GC content, and adapter contamination.
     - Users can review these files to ensure the input reads meet quality standards before downstream analysis.

### 3. **Variant Calling Results**
   - **Location**: `results/variant_calls/`
   - **File**: `final.vcf`
   - **Description**:
     - This VCF file contains all raw variants detected during the variant calling step.
     - Includes potential single nucleotide polymorphisms (SNPs), insertions, and deletions (indels).
     - May contain some low-confidence variants, which are filtered in the `filtered.vcf` file for more accurate downstream analyses.
    
### 4. **Filtered VCF Files**
   - **Location**: `results/filtered_vcf/`
   - **File**: `filtered.vcf`
   - **Description**:
     - This Variant Call Format (VCF) file contains high-confidence variants that meet the pipeline's filtering criteria.
     - Common filters applied include minimum read depth and quality thresholds, ensuring only reliable variants are included.
     - Used for downstream analysis such as annotation or population-level comparisons.


---

### Temporary Output in Work Directory

These intermediate files are stored in the `work/` directory and used during the pipeline execution:

1. **Genome Index Files**
   - **Location**: Intermediate files in the `work/` directory.
   - **Description**: Includes SAMtools index (`*.fai`), Picard sequence dictionary (`*.dict`), and STAR genome index files.

2. **BAM Files**
   - **Location**: Intermediate files in the `work/` directory.
   - **Description**: Aligned reads saved in `Aligned.sortedByCoord.out.bam`.

3. **Split BAM File**
   - **Location**: `split.bam` in the `work/` directory.
   - **Description**: Handles intron-spanning alignments for RNA-Seq-specific challenges.

## Expected Output Directory Tree
![image](https://github.com/user-attachments/assets/e0fd635d-4476-4904-85e0-80bb511d5d9b)

•	data/: Contains all input files.

•	nextflow.config: Configuration file specifying Docker containers and environment settings.

•	results/: Stores output files, organized into subdirectories for clarity.

•	work/: Contains temporary files created during execution.

•	variant_calling.nf: Main Nextflow pipeline script orchestrating the entire analysis.

---

## Conclusion 

The Variant Calling Analysis Pipeline offers a comprehensive, efficient, and reproducible solution for processing RNA-Seq data, from raw reads to high-confidence variant detection. By integrating industry-standard tools such as FastQC, SAMtools, Picard, STAR, and GATK within a Dockerized Nextflow environment, the pipeline ensures consistency and eliminates software compatibility issues. Its automated workflow simplifies complex bioinformatics tasks, such as quality control, alignment, variant calling, and filtering, while providing outputs that are ready for downstream analyses like annotation or comparative studies. Whether run locally or in a pre-configured Gitpod environment, this pipeline is well-suited for both large-scale and smaller projects requiring reliable and reproducible RNA-Seq data analysis.

---

## References 
[1] Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Software]. Available at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[2] Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

[3] Broad Institute. (2019). Picard Toolkit [Software]. Broad Institute. GitHub Repository. Available at: https://broadinstitute.github.io/picard/

[4] Dobin, A., Davis, C. A., Schlesinger, F., et al. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635

[5] McKenna, A., Hanna, M., Banks, E., et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

[6] Van der Auwera, G. A., & O’Connor, B. D. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra. O’Reilly Media.

## Contact

For questions or issues, please contact:
	•	Name: Yiyang Wang
	•	Email: yiyang.wang2002@gmail.com
	•	GitHub: yiyang2002 write one paragraph of conclusion for this pipieline
