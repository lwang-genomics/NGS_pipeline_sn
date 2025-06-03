## II. ChIP-Seq Processing Pipeline 

This Snakemake pipeline streamlines the analysis of ChIP-seq data by automating key processing steps for both single-end and paired-end sequencing. Designed for robustness and flexibility, it handles raw FASTQ files through quality control, optional adapter trimming, genome alignment, peak calling, and signal track generation. With a modular structure and centralized configuration via config.yaml, it is well-suited for reproducible, multi-sample ChIP-seq projects across various experimental designs.

### Features

- Automatic detection of input FASTQ files (\*.R1.fq.gz, \*.R1.fastq, etc.)
- Supports both paired-end and single-end reads
- Optional read trimming with **Trimmomatic**
- Read alignment with **BWA**
- Filtering, sorting, and indexing via **SAMtools**
- Peak calling with **MACS2** (configurable for narrow or broad peaks)
- Generation of normalized BigWig files with **deepTools**
- Integrated summary with **MultiQC**
- Possible to provide a customized species genome
- Configurable via a single config.yaml

### Requirements
	Snakemake
	BWA
	SAMtools
	Trimmomatic
	FastQC
	deepTools
	MACS2
	MultiQC

Install dependencies using mamba, conda, or your preferred environment manager.

### Usage
1. Copy and paste chip_seq.smk and config.yaml into your working folder with all sample FASTQ files


2. Modify the config.yaml file as needed


3. Run Snakemake 
```
snakemake -s atac_seq.smk --cores 8
```
Optional: Log output to a file and run in background:

```
snakemake -s atac_seq.smk --cores 8 > snakemake.log 2>&1 &
```


### Configuration

Edit the config.yaml file to customize:

```text
threads: 4                       # Threads per rule
mapq: 5                          # Minimum MAPQ for filtering
read_type: paired                # Read type: paired or single
genome: 
	name: hg38                   # Genome name (e.g., hg38 or mm10)
	bwa_index:[bwa_index/genome] # Path to BWA index (+ prefix)
peaktype: narrow                 # Peak type: narrow or broad
keep_intermediate: false         # Whether to retain intermediate files
skip_trimming: false             # Skip trimming step if reads are already preprocessed
```


### Output Files
- sampleX_filtered_sorted.bam(.bai) – Filtered, sorted, and indexed BAM files
- sampleX.bw – BigWig signal track normalized by CPM
- sampleX_peaks.narrowPeak or .broadPeak – Peak calls from MACS2
- multiqc_report.html – Combined QC report


### Example Project Structure

```
project/
├── chip_seq.smk
├── config.yaml
├── sample1.R1.fq.gz
├── sample1.R2.fq.gz
├── sample2.R1.fq.gz
├── sample2.R2.fq.gz
└── ...
```




## III. ATAC-Seq Processing Pipeline 

This Snakemake pipeline provides a lightweight, scalable, and reproducible solution for processing ATAC-seq data. It automatically detects paired-end FASTQ files with flexible naming patterns and processes them through trimming, alignment, mitochondrial read removal, filtering, peak calling, QC, and visualization. Configurable via a single config.yaml file, it is ideal for streamlined, multi-sample ATAC-seq analysis in any computing environment.


### Features

- Automatic detection of input FASTQ files (*.R1.fq.gz, *.R1.fastq, etc.)
- Configurable trimming step (can be toggled on/off)
- Read alignment with **BWA**
- Mitochondrial reads removal, filtering and sorting via **SAMtools**
- Peak calling with **MACS2** (configurable for narrow or broad peaks)
- Generation of normalized BigWig files with **deepTools**
- ATAC-specific QC with **ATAQV**
- Integrated summary with **MultiQC**
- Possible to provide a customized species genome
- Configurable via a single config.yaml

### Requirements
	Snakemake
	BWA
	SAMtools
	Trimmomatic
	FastQC
	deepTools
	MACS2
	ataqv
	MultiQC

Install dependencies using mamba, conda, or your preferred environment manager.

### Usage
1. Copy and paste atac_seq.smk and config.yaml into the folder that stores all sample fastq files


2. Modify the config.yaml file as needed


3. Run snakemake 
```
snakemake -s atac_seq.smk --cores 8
```
Optional: Log output to a file and run in background:

```
snakemake -s atac_seq.smk --cores 8 > snakemake.log 2>&1 &
```


### Configuration

Edit the config.yaml file to customize:

```text
threads: 4                   # Threads per rule
mapq: 5                      # Minimum MAPQ for filtering
genome: hg38                 # Genome name (e.g., hg38 or mm10)
bwa_index:[bwa_index/genome] # Path to BWA index (+ prefix)
peaktype: narrow             # Peak type: narrow or broad
keep_intermediate: false     # Whether to retain intermediate files
skip_trimming: false         # Skip trimming step if reads are already preprocessed
```


### Output Files
- sampleX_filtered_sorted.bam(.bai) – filtered, sorted, and indexed BAM files
- sampleX.bw – normalized BigWig files
- sampleX_peaks.narrowPeak or broadPeak – peak calls
- sampleX_ataqv_metrics.json – QC metrics
- multiqc_report.html – summary of QC and analysis


### Example Project Structure

```
project/
├── atac_seq.smk
├── config.yaml
├── sample1.R1.fq.gz
├── sample1.R2.fq.gz
├── sample2.R1.fq.gz
├── sample2.R2.fq.gz
└── ...
```

### License

MIT License

### Acknowledgments

This pipeline integrates many excellent open-source bioinformatics tools. Credit goes to the developers of Snakemake, BWA, SAMtools, MACS2, deepTools, MultiQC and so on.



Please let me know if you have any questions or suggestions about my pipeline tool!



