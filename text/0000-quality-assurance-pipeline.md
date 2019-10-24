# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
  - [Current Process](#current-process)
  - [Important Metrics](#important-metrics)
  - [Proposals](#proposals)
- [Specification](#specification)
  - [Workflow Description](#workflow-description)
- [Items Still In-Progress](#items-still-in-progress)
- [Outstanding Questions](#outstanding-questions)

# Introduction

This RFC documents an automated pipeline workflow for vetting St. Jude Cloud genomic data, covering both existing data and new data uploads to the platform. The end goal is to publish results from various tools, but currently we hope to discuss which quality metrics and statistics are important to the community as a whole.

# Motivation

Since introducing Real-Time Clinical Genomics, there is a need for an automated quality assurance pipeline that guarantees uploaded data meets predefined standards.  Guaranteeing the data integrity and the reproducibility of these results allows St. Jude to publish statistics that are of interest to scientists and researchers.

The ultimate goal is to present a comprehensive report much like the [example MultiQC report](https://multiqc.info/examples/rna-seq/multiqc_report.html) for each dataset and sequencing type (and ideally, also on a sample level). This would make the quality of data offered to researchers and scientists accessible. We hope this RFC becomes a forum for open community discussion of quality properties and attributes that are helpful and practical.

# Discussion

The quality metrics discussed here are sequence and mapping quality metrics.   Other metrics related to nucleic acid integrity or library quality are not part of the current process.  

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

- Per Base Sequence Quality ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html))

The "Per Base Sequence Quality" module from FastQC shows the distribution of quality scores across all bases at each position in the reads. In our case, this is just to inform our end users — the quality of the sequencing run has already been assessed by the lab upstream, so there is no changing it at this point.

- Overrepresented Sequences ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html))

The "Overrepresented Sequences" module from FastQC displays sequences (at least 20bp) that occur in more than 0.1% of the total number of sequences and will help identify contamination (vector, adapter sequences, etc.).

- Reads Genomic Origin ([Qualimap](http://qualimap.bioinfo.cipf.es/))

The "Reads Genomic Origin" from Qualimap determines how many alignments fall into exonic, intronic, and intergenic regions. Even if there is a high genomic mapping rate, it is necessary to check where the reads are being mapped. It should be verified that the mapping to intronic regions and exons are within acceptable ranges. Abnormal results could indicate issues such as DNA contamination.

- rRNA Content (?)

Verify that excess ribosomal content is filtered/normalized across samples to ensure that alignment rates and subsequent normalization of data is not skewed.

- Transcript Coverage and 5’-3’ Bias ([Qualimap](http://qualimap.bioinfo.cipf.es/))

Libraries prepared with polyA selection may have higher biased expression in 3’ region. If reads primarily accumulate at the 3’ end of transcripts (in poly(A)-selected samples), this might indicate the starting RNA was of low quality.

- Junction Analysis ([Qualimap](http://qualimap.bioinfo.cipf.es/))

Analysis of known, partly known, and novel junction positions in spliced alignments.

- Strand Specificity ([RSeQC](http://rseqc.sourceforge.net/))

Verification/sanity check of how reads were stranded for the RNA sequencing (stranded or unstranded protocol). 

- GC Content Bias (?)

GC profiles are typically remarkably stable. Even small/minor deviations could indicate a problem with the library used (or bacterial contamination).

## Proposals

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
| `fastqc`                 | Screen for GC content and adapter contamination            |
| `qualimap bamqc`         | Screen for mapping quality, coverage, and duplication rate |
| `qualimap rnaseq`        | Screen for RNA-Seq bias and junction analysis              |
| `rseqc infer_experiment` | Determine RNA-SEQ strandedness and reads                   |
| `multiqc`                | Report aggregation                                         |

Note: Specific options such as memory size thresholds and thread count have been left out.

# Items Still In-Progress

- [ ] Analysis tools for other types of sequencing (ChIP seq)
- [ ] Useful metadata from various stages (sample collection, laboratory, pre-sequencing, sequencing, post-sequencing)

# Outstanding Questions

- What thresholds or metrics differentiate a poor-quality sample from a high-quality one?
- What other metrics or properties would be valuable?
- What is best way to define and handle outliers?
- What is the best way to examine cohort integrity, meaning category-based tests of samples to find experimental outliers that are of sufficient quality if examined alone? 
