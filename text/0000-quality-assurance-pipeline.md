# Table of Contents <!-- omit in toc -->

- [Introduction](#Introduction)
- [Motivation](#Motivation)
- [Current Process](#Current-Process)
- [Proposals](#Proposals)
- [Workflow Description](#Workflow-Description)
- [Items Still In-Progress](#Items-Still-In-Progress)
- [Outstanding Questions](#Outstanding-Questions)

# Introduction

This RFC seeks to establish an automated pipeline workflow around how genomic data on St. Jude Cloud is vetted, covering both existing data and new uploads to the platform. The end goal for this would be to publish results from various tools, but it hopes to draw discussion around what metrics and statistics are important to the community as a whole.

# Motivation

With the introduction of Real-Time Clinical Genomics, there exists a need for an automated quality assurance pipeline guaranteeing any uploaded data meets predefined standards. By guaranteeing the integrity of our data and the reproducibility of these results, it would allow St. Jude to publish statistics about the genomics data hosted on our platform that might be of interest to other scientists and researchers.

The ultimate goal is to be able to present a comprehensive report much like the [example MultiQC report](https://multiqc.info/examples/rna-seq/multiqc_report.html) separated by dataset and sequencing type (and ideally, also on a sample level). This would aid in visibility into the quality and type of data hosted to researchers and scientist. This RFC hopes to present a form for open discussion to the community regarding what type of other properties/attributes would be helpful and practical.

# Current Process

Because St. Jude Cloud currently provides three-platform whole-genome (WGS), whole-exome (WES), and transcriptome (RNA-Seq) sequencing data, it is important to differentiate how we currently run our current quality control workflow on each type of sequencing.

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

# Important Metrics

- Per Base Sequence Quality ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html))

The "Per Base Sequence Quality" module from FastQC will show the distribution of quality scores across all bases at each position in the reads. It will automatically determine the encoding method used, but this should be cross-referenced with the actual encoding method.

- Overrepresented Sequences ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html))

The "Overrepresented Sequences" module from FastQC displays sequences (at least 20bp) that occur in more than 0.1% of the total number of sequences and will help identify any sort of contamination (vector, adapter sequences, etc.).

- Reads Genomic Origin ([Qualimap](http://qualimap.bioinfo.cipf.es/))

The "Reads Genomic Origin" from Qualimap is able to determine how many alignments fall into exonic, intronic, and intergenic regions. Even if there is a high genomic mapping rate, it is necessary to check where the reads are being mapped to. It should be verified that the mapping to intronic regions and exons are within acceptable ranges. Any abnormal results could indicate issues such as DNA contamination.

- rRNA Content (?)

Verify that excess ribosomal content is filtered/normalized across samples to ensure that alignment rates and subsequent normalization of data is not skewed.

- Transcript Coverage and 5’-3’ Bias ([Qualimap](http://qualimap.bioinfo.cipf.es/))

Libraries prepared with polyA selection have the possibility to lead to high expression in 3’ region. If reads primarily accumulate at the 3’ end of transcripts (in poly(A)-selected samples), this might indicate the starting material was of low RNA quality.

- Junction Analysis ([Qualimap](http://qualimap.bioinfo.cipf.es/))

Analysis of known, partly known, and novel junction positions in spliced alignments.

- Strand Specificity ([RSeQC](http://rseqc.sourceforge.net/))

Verification/sanity check of how reads were stranded for the RNA sequencing (stranded or unstranded protocol). 

- GC Content Bias (?)

GC profiles are typically remarkably stable. Even small/minor deviations could indicate a problem with the library used (or bacterial contamination).

# Proposals

- Add [`RSeQC v3.0.0`](http://rseqc.sourceforge.net), specifically [`infer_experiment`] and [`junction_annotation`].

[`infer_experiment`]: http://rseqc.sourceforge.net/#infer-experiment-py
[`junction_annotation`]: http://rseqc.sourceforge.net/#junction-annotation-py

- Include md5 hash as an annotation property for vended files.

# Workflow Description

The end workflow (covering both our current process and the addition of the new tool) would be as following:

| Command                     | Purpose                                                    |
| --------------------------- | ---------------------------------------------------------- |
| `samtools quickcheck`       | Validate BAM headers and EOF block existence               |
| `md5sum`                    | For comparison to md5 vended file property                 |
| `picard ValidateSamFile`    | Ensure validity of file                                    |
| `samtools flagstat`         | Generate flag statistics                                   |
| `fastqc`                    | Screen for GC content and adapter contamination            |
| `qualimap bamqc`            | Screen for mapping quality, coverage, and duplication rate |
| `qualimap rnaseq`           | Screen for RNA-Seq bias and junction analysis              |
| `rseqc infer_experiment`    | Determine RNA-SEQ strandedness and reads                   |
| `rseqc junction_annotation` | Compare detected splice junctions to reference gene model  |
| `multiqc`                   | Report aggregation                                         |

Note: Specific options such as memory size thresholds and thread count have been left out.

# Items Still In-Progress

- [ ] Analysis tools for other types of sequencing
- [ ] Useful metadata from various stages (sample collection, laboratory, pre-sequencing, sequencing, post-sequencing)

# Outstanding Questions

- What thresholds or metrics differentiate a poor-quality sample from a high-quality one?
- What other metrics or properties would be valuable?
