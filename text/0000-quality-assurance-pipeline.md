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
