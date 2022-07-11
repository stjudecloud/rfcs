# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)

# Introduction

This RFC lays out the specification for the scRNA-Seq mapping pipeline.

# Motivation

To provide the  community access to data from scRNA-Seq experiments performed at St. Jude Children's Research Hospital, we propose the following data harmonization pipeline. The goal of this pipeline is to provide harmonized alignments for scRNA-Seq data. For this pipeline, we will make no recommendations on downstream analysis, focusing instead on harmonizing the underlying sequencing data and leaving analysis decisions to the user.

# Discussion

The development of bulk RNA-Seq enabled efficient profiling of transcripts in a sample. The ability to interrogate the transcriptome in an unbiased manner, largely replaced previous methods such as microarrays or RT-qPCR. A bulk RNA-Seq experiment sample is a mixture of cells. This method is useful for analyses such as comparing expression between tumor and normal samples. However, this limits the analysis to the average expression over a population of cells and is unable to capture the heterogeneity of gene expression across individual cells. Newer methods have enabled RNA-Seq to be applied to single cells. This enables new questions about cell-specific transcriptome changes. This methodology has been employed in large projects such as the [Human Cell Atlas](https://www.humancellatlas.org/) to capture the cellular diversity within organisms.

## scRNA-Seq methodology

Generally, scRNA-Seq differs somewhat from bulk RNA-Seq. Most approaches generate three key pieces of information: 1. cDNA fragment from the RNA transcript, 2. a cell barcode that identifies in which cell the RNA was expressed, and 3. a unique molecular identifier (UMI) to enable collapse of PCR duplicates.

There are several commerical options for performing scRNA-Seq. One such option is the 10x Chromium platform. In their approach using paired-end reads, the celluar barcode and UMI are found in read 1 and the transcript sequence is found in read 2.

The canonical tool for processing 10x Chromium scRNA-Seq data is their [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) tool. This tool uses STAR to perform alignment. It then uses the transcriptome annotation GTF to label reads as exonic, intronic, or intergenic. It then performs some additional refinement on the alignments to determine which reads are transcriptomic. These reads are then carried into UMI counting. The UMI counting step produces a raw gene by cell matrix of counts.

## Processing Pipeline

For processing of bulk RNA-Seq data, a common standard workflow is to align against the genome using [STAR](https://github.com/alexdobin/STAR). Quantification is then done using a package such as [HTSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html). An example of a bulk RNA-Seq pipeline is the [St. Jude Cloud RNA-Seq standard workflow v2.0.0](0001-rnaseq-workflow-v2.0.0.md).

### Aligner choice

For consistency with other St. Jude Cloud datasets, we are not considering psuedoalignment methods such as [kallisto](https://www.kallistobus.tools/) and [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html). This leaves the previously described Cell Ranger method and STARsolo as the primary candidates for use in St. Jude Cloud.

STARsolo is a standalone pipeline that is integrated into the STAR RNA-Seq aligner. It was developed with a goal to generate highly similar results to those of Cell Ranger. STARsolo has the advantage of typically being significantly faster than Cell Ranger on the same dataset. It reimplments the algorithms of Cell Ranger in order to produce highly similar results. [Bruning, et. al.](https://doi.org/10.1093/gigascience/giac001) (2022) provide a useful comparison of commonly used scRNA-Seq methods.

<img src="../resources/0002-scRNA-Seq-workflow/aligner_comparison.jpeg" width=800">

### Quantification choice

Quantification is a product of both the Cell Ranger and STARsolo packages. Quantification could also be accomplished using a pseudoalignment method, such as Kallisto or Alevin, if there was no desire to share the raw data in BAM format. However that would be a departure from the normal St. Jude Cloud methodology, so that will not be discussed here. Therefore, there is no additional quantification choice to be made beyond that of the aligner package.

Using either STARsolo or Cell Ranger for analysis requires selection of reference files. For the genome, we will continue using GRCh38_no_alt. The transcript annotation GTF is more complicated. 10x [provides](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) several versions of the annotation files used in Cell Ranger. However for consistency, we would prefer to use GENCODE v31 for consistency with our other annotation-based methods.