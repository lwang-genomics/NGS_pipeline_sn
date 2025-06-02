## II. ATAC-Seq Processing Pipeline 

This repository contains a modular and automated Snakemake pipeline for processing ATAC-seq data. Designed to be scalable, reproducible, and easy to configure, the pipeline takes raw FASTQ files and produces filtered BAMs, normalized BigWig tracks, quality metrics, and peak calls in a streamlined manner. It supports configuration through a config.yaml file, with flexible options for trimming, genome setup, peak type, and file retention.

### Features

- Automatic detection of input FASTQ files (supports .fastq, .fq, and .gz)
- Optional trimming using Trimmomatic
- Quality check using FastQC
- Read alignment with BWA
- Filtering and sorting via SAMtools
- Optional removal of mitochondrial reads
- Peak calling with MACS2, supporting both narrow and broad peaks
- Generation of normalized BigWig files with deepTools
- Quality control summary with ataqv and MultiQC
- Temporary intermediate file cleanup using Snakemake’s temp()
- Configurable via a single config.yaml

### Quick Start

1. Clone the repo
```
git clone https://github.com/yourname/NGS_pipeline_sn.git
cd NGS_pipeline_sn

# Prepare your input directory with FASTQ files (e.g., sample1.R1.fq.gz, sample1.R2.fq.gz, ...)
# Customize parameters in config.yaml as needed

# Run the pipeline (from within the FASTQ directory or adjust FASTQ_DIR in the script)
snakemake --cores 8
```

### Required Configuration (config.yaml)

```text
threads: 4                 # Threads per rule
mapq: 5                    # Minimum MAPQ for filtering
genome: hg38               # Genome name (e.g., hg38 or mm10)
genome_size: hs            # MACS2 genome size flag (e.g., hs for human, mm for mouse)
bwa_index: ref/genome      # Path to BWA index (prefix only)
peaktype: narrow           # Peak type: narrow or broad
keep_intermediate: false   # Whether to retain intermediate files
skip_trimming: false       # Skip trimming step if reads are already preprocessed
```


### Output Files
- sample_filtered_sorted.bam(.bai) – filtered, sorted, and indexed BAM files
- sample.bw – normalized BigWig files
- sample_peaks.narrowPeak or broadPeak – peak calls
- sample_ataqv_metrics.json – QC metrics
- multiqc_report.html – summary of QC and analysis

### Requirements
	•Snakemake
	•BWA
	•SAMtools
	•Trimmomatic
	•FastQC
	•deepTools
	•MACS2
	•ataqv
	•MultiQC

### Notes
	•Designed to be run from within the folder containing input FASTQ files.
	•Output file names are sample-name prefixed and automatically inferred.
	•Intermediate files can be retained by setting keep_intermediate: true.




