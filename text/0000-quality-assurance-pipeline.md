# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
  - [Current Process](#current-process)
  - [Important Metrics](#important-metrics)
  - [Thresholds and Metrics for Specific Applications](#thresholds-and-metrics-for-specific-applications)
    - [Metrics for WES](#metrics-for-wes)
    - [Metrics for RNAseq](#metrics-for-rnaseq)
    - [Proposals](#proposals)
- [Specification](#specification)
  - [Workflow Description](#workflow-description)
- [The St Jude Genomics QC Process](#the-st-jude-genomics-qc-process)
  - [Installation](#installation)
      - [Anaconda Environment](#anaconda-environment)
    - [Samtools Quickcheck](#samtools-quickcheck)
    - [Picard ValidateSamFile](#picard-validatesamfile)
    - [Sambamba Flagstat](#sambamba-flagstat)
    - [FASTQC](#fastqc)
    - [Qualimap Bam QC](#qualimap-bam-qc)
    - [Qualimap RNA seq QC](#qualimap-rna-seq-qc)
    - [md5sum](#md5sum)
    - [RseQC infer_experiment.py](#rseqc-infer_experimentpy)
    - [Report Aggregation](#report-aggregation)
- [Items Still In-Progress](#items-still-in-progress)
- [Outstanding Questions](#outstanding-questions)

# Introduction

This RFC documents an automated pipeline workflow for vetting St. Jude Cloud genomic data, covering both existing and new data uploads to the platform.  The end goal is to publish the quality control results from various tools.  But currently, we hope to discuss which quality metrics and statistics are important to the bioinformatics community.  Further, we invite the community to comment on the best methods for selecting thresholds for quality metrics. 

# Motivation

Since introducing Real-Time Clinical Genomics, we need an automated quality assurance pipeline that guarantees uploaded data meets predefined standards.   Guaranteeing data integrity and the reproducibility of these results allows St. Jude to assure scientists and researchers that the data we provide is useful.

The ultimate goal is to present a comprehensive report much like the example [ MultiQC report](https://multiqc.info/examples/rna-seq/multiqc_report.html) for each dataset and sequencing type (and ideally, also on a sample level).  This would make the quality of data offered to researchers and scientists accessible.  We hope this RFC becomes a forum for open community discussion of quality properties and attributes that are helpful and practical.

# Discussion

The quality metrics discussed here are sequence and mapping quality metrics.   Other metrics related to nucleic acid integrity or library quality are not part of the current process.  Pre-sequencing quality metrics, however, are clearly important and part of our long term interests.  Thus, we invite comments on those metrics as well.


## Current Process

Currently, St. Jude Cloud provides three sequencing data types: whole-genome (WGS), whole-exome (WES), and transcriptome (RNA-Seq) data.  It is important to differentiate our quality control workflows for each type of sequencing.

Our current process to vet and screen data consists of the following tools:

| Tool                     | Version   |
| ------------------------ | --------- |
| `samtools flagstat`      | [v1.9]    |
| `fastqc`                 | [v0.11.8] |
| `qualimap bamqc`         | [v2.2.2]  |
| `qualimap rnaseq`        | [v2.2.2]  |
| `picard ValidateSamFile` | [v2.20.2] |
| `multiqc`                | [v1.7]    |

[v1.9]: http://www.htslib.org/doc/samtools.html
[v0.11.8]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[v2.2.2]: http://qualimap.bioinfo.cipf.es/doc_html/command_line.html
[v2.20.2]: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_ValidateSamFile.php
[v1.7]: https://multiqc.info/

## Important Metrics

- Percent Aligned  

Also known as mapping percentage, this indicator of quality, when high, verifies the mapping process/genome was correct and is consisitent with sample purity.

- Per Base Sequence Quality ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html))

The "Per Base Sequence Quality" module from FastQC shows the distribution of quality scores across all bases at each position in the reads.  In our case, this is just to inform our end users — the quality of the sequencing run has already been assessed by the lab upstream.  So, there is no changing it at this point.

- Overrepresented Sequences ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html))

The "Overrepresented Sequences" module from FastQC displays sequences (at least 20bp) that occur in more than 0.1% of the total number of sequences and will help identify contamination (vector, adapter sequences, etc.).

- Reads Genomic Origin ([Qualimap](http://qualimap.bioinfo.cipf.es/))

The "Reads Genomic Origin" from Qualimap determines how many alignments fall into exonic, intronic, and intergenic regions.  Even if there is a high genomic mapping rate, it is necessary to check where the reads are being mapped.  It should be verified that the mapping to intronic regions and exons are within acceptable ranges.  Abnormal results could indicate issues such as DNA contamination.

- rRNA Content (?)

Verify that excess ribosomal content is filtered/normalized across samples to ensure that alignment rates and subsequent normalization of data is not skewed.

- Transcript Coverage and 5’-3’ Bias ([Qualimap](http://qualimap.bioinfo.cipf.es/))

Libraries prepared with polyA selection may have higher biased expression in 3’ region.  If reads primarily accumulate at the 3’ end of transcripts (in poly(A)-selected samples), this might indicate the starting RNA was of low quality.

- Junction Analysis ([Qualimap](http://qualimap.bioinfo.cipf.es/))

Analysis of known, partly known, and novel junction positions in spliced alignments.

- Strand Specificity ([RSeQC](http://rseqc.sourceforge.net/))

Verification/sanity check of how reads were stranded for the RNA sequencing (stranded or unstranded protocol). 

- GC Content Bias (?)

GC profiles are typically remarkably stable.  Even small/minor deviations could indicate a problem with the library used (or bacterial contamination).

## Thresholds and Metrics for Specific Applications 

 To apply quality control metrics to vett data, we need reasonable thresholds that are practically acheivable and neither too lax or too strict.  Our preference is for statistically or empirically determined thresholds rather than arbitrary estimates.  By statistical thrresholds, we are referring to distributional tests that formally define outliers.  By empirical thresholds, we are referring to standards below which data analysis or interpretation are degraded.   Statistical tests can be performed on large populations of QC data. We are already in postion to do that today.  Empirical tests, however, require foreknowledge of the correct results.  This requires experimental design and implementation through a laboratory at some cost.## Metrics for WGS
    
The quality metrics of special concern for WGS include depth of coverage and genomic regional coverage. Mapping quality is also critical.  The analysis of whole genome sequencing to call variants depends on depth and sample purity. Accurate calls are made through replication and contamination creates false positives.  So metrics that are sensitive to impurity are valuable.

### Metrics for WES

The quality metrics of special concern for WES include depth of coverage in exomic regions. Mapping quality, % mapped and duplication rate are also important.

### Metrics for RNAseq

The quality metrics of special concern for RNAseq include mapping percentage, percentage properly paired reads, and exomic regional coverage.  Mapping quality is also critical.


### Proposals

- Add [`RSeQC v3.0.0`](http://rseqc.sourceforge.net), specifically [`infer_experiment`].

[`infer_experiment`]: http://rseqc.sourceforge.net/#infer-experiment-py

- Include md5 hash as an annotation property for vended files.

# Specification

## Workflow Description

The end workflow (covering both our current process and the addition of the new tool) would be as following:

| Command                  | Purpose                                                    |
| ------------------------ | ---------------------------------------------------------- |
| `samtools quickcheck`    | Validate BAM headers and EOF block existence               |
| `md5sum`                 | For comparison to md5 vended file property                 |
| `picard ValidateSamFile` | Ensure validity of file                                    |
| `samtools flagstat`      | Generate flag statistics                                   |
| `sambamba flagstat`      | Generate flag statistics using a samtools alternative      |
| `fastqc`                 | Screen for GC content and adapter contamination            |
| `qualimap bamqc`         | Screen for mapping quality, coverage, and duplication rate |
| `qualimap rnaseq`        | Screen for RNA-Seq bias and junction analysis              |
| `rseqc infer_experiment` | Determine RNA-SEQ strandedness and reads                   |
| `multiqc`                | Report aggregation                                         |

Note: Specific options such as memory size thresholds and thread count are below.

# The St Jude Genomics QC Process

These are generic instructions for running each of the tools in our pipeline.  We run our pipeline in a series of QC scripts that are tailored for our compute cluster, so those commands may not apply elsewhere.  Instead we've supplied examples of the commands used to each package.  Our default memory is 80G and we employ 4 threads for these processes.

## Installation

We presume anaconda is available and installed. If not please follow the link to [anaconda](https://www.anaconda.com/) first.

#### Anaconda Environment


```bash
conda create --name bio-qc \
    --channel bioconda \
    fastqc==0.11.8 \
    picard==2.20.2  \
    qualimap==2.2.2c \
    rseqc==3.0.0  \
    sambamba==0.6.6 \
    samtools==1.9  \
    -y

conda activate bio-qc
```
### Samtools Quickcheck

Very basic BAM file validation:

```bash
# Should not be relied on for file integrity, only checks for header and EOF
samtools quickcheck -v *.bam > bad_bams.txt && echo "all ok" || echo "some files failed check, see bad_bams.txt"
```
### Picard ValidateSamFile

BAM file validation:

```bash
# This method is used to assess file integrity

picard ValidateSamFile \
    I=$BAM \                                        # specify bam file
    MODE=SUMMARY\                                   # concise output
    INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE \   # lower stringency faster processing time
    OUTPUT=$OUTDIR/$BAM_BN.validate.txt             # output directory and file
```

### Sambamba Flagstat

Summary statistics of read counts and mapping status

```bash
# Sambamba includes faster reliable implementation of samtools commands  

sambamba flagstat -t $NUM_THREADS $BAM \ # number of threads and bam filename
         > $OUTDIR/$BAM_BN.flagstat.txt  # output directory and file
```

### FASTQC

Standard sequence quality check:

```bash

fastqc $BAM \ #  bam filename
     -o $OUTDIR   # output directory
        
```
 
### Qualimap Bam QC

Comprehensive QC statistics includes read stats, coverage, mapping quality, insert size, mismatches etc.

```bash

qualimap bamqc -bam $BAM \          # bam filename
    --java-mem-size=$MEM_SIZE \     # memory 
    -nt $NUM_THREADS \              # threads requested
    -nw 400 \                       # number of windows
    -outdir $QBAMQC_OUT             # output directory
```
### Qualimap RNA seq QC

Comprehensive QC statistics tailored for RNA seq files.

```bash

qualimap rnaseq -bam $BAM \         # bam filename
    -gtf $GTF_REF                   # transcript definition file
    --java-mem-size=$MEM_SIZE \     # memory 
    -pe                             # specify paired end
    -outdir $QBAMQC_OUT             # output directory
```

### md5sum

Check size and integrity of files.

```bash

md5sum $BAM \                       # bam filename
    > $OUTDIR/$BAM.md5              # output directory
           
```

### RseQC infer_experiment.py

Python script that tests for strandedness.

```bash

infer_experiment.py -i $BAM \       # bam filename
                -r $BED_REF \       # reference in bed format
                > $OUTDIR/$BAM_BN.infer_experiment.txt # output directory and filename
           
```

### Report Aggregation

Package combines output from other QC tools in easily reviewed html format.

```bash
multiqc /path/to/outdir
```


# Items Still In-Progress

- [ ] Analysis tools for other types of sequencing (ChIP seq)
- [ ] Useful metadata from various stages (sample collection, laboratory, pre-sequencing, sequencing, post-sequencing)

# Outstanding Questions

- What thresholds or metrics differentiate a poor-quality sample from a high-quality one?
- What other metrics or properties would be valuable?
- What is best way to define and handle outliers?
- What is the best way to examine cohort integrity? This means experimental category-based tests of samples to find outliers that are of sufficient quality if examined alone.  Outliers in this case may indicate classification errors or rare biological conditions.  Which metrics are best tested here?
