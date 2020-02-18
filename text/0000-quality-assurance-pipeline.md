# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
- [Specification](#specification)
- [Items Still In-Progress](#items-still-in-progress)
- [Outstanding Questions](#outstanding-questions)

# Introduction

This RFC documents an automated workflow for assessing the integrity and quality
of St. Jude Cloud genomics data. The end goal is to publish a collection of
metrics that users can leverage to assess the quality of the data available.
Furthermore, we outline the method used internally to vet the data before we
publish it. You can find the relevant discussion on the [associated pull request](https://github.com/stjudecloud/rfcs/pull/3).

# Motivation

Since the introduction of uploading clinical genomics data in real-time (the
"Real-Time Clinical Genomics" initiative), we need an automated quality
assurance pipeline that guarantees uploaded data meets predefined standards.
Guaranteeing data integrity and the reproducibility of these results allows St.
Jude to assure scientists and researchers that the data we provide is useful. As
much as possible, we'd like to automate this process to ensure it scales. The end-goal is to present a comprehensive report, much like the example [
MultiQC report](https://multiqc.info/examples/rna-seq/multiqc_report.html), for
each dataset + sequencing type tuple.

# Discussion

The quality metrics discussed here are sequence and mapping quality metrics.
Other metrics related to nucleic acid integrity or library quality are typically
done upstream in the genomics lab contributing the data. Pre-sequencing quality metrics, however, are clearly important and part of our long term interests.

## Current Process

Currently, St. Jude Cloud provides three sequencing data types: whole-genome
(WGS), whole-exome (WES), and transcriptome (RNA-seq) data.  It is important to
differentiate our quality control workflows for each type of sequencing. At the
time that the RFC was authored, we use the following tools to _manually_ vet
data.

> Note that one of the goals of this RFC is to constantly update these
versions, so they are likely not what will be used in the pipeline.

| Tool                     | Version | Website          |
| ------------------------ | ------- | ---------------- |
| `samtools flagstat`      | v1.9    | [link](samtools) |
| `fastqc`                 | v0.11.8 | [link](fastqc)   |
| `qualimap bamqc`         | v2.2.2  | [link](qualimap) |
| `qualimap rnaseq`        | v2.2.2  | [link](qualimap) |
| `picard ValidateSamFile` | v2.20.2 | [link](picard)   |
| `multiqc`                | v1.7    | [link](multiqc)  |


## Important Metrics

| Name                               | Produced By          | Description                                                                                                                                                                                                                                                                                                                                                                                                  |
| ---------------------------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| % Aligned                          | [Samtools][samtools] | Also known as mapping percentage, this indicator of quality, when high, verifies the mapping process/genome was correct and is consisitent with sample purity.                                                                                                                                                                                                                                               |
| Per Base Sequence Quality          | [FastQC][fastqc]     | The "Per Base Sequence Quality" module from FastQC shows the distribution of quality scores across all bases at each position in the reads.  In our case, this is just to inform our end users — the quality of the sequencing run has already been assessed by the lab upstream.  So, there is no changing it at this point.                                                                                |
| Overrepresented Sequences          | [FastQC][fastqc]     | The "Overrepresented Sequences" module from FastQC displays sequences (at least 20bp) that occur in more than 0.1% of the total number of sequences and will help identify contamination (vector, adapter sequences, etc.).                                                                                                                                                                                  |
| Reads Genomic Origin               | [Qualimap][qualimap] | The "Reads Genomic Origin" from Qualimap determines how many alignments fall into exonic, intronic, and intergenic regions.  Even if there is a high genomic mapping rate, it is necessary to check where the reads are being mapped.  It should be verified that the mapping to intronic regions and exons are within acceptable ranges.  Abnormal results could indicate issues such as DNA contamination. |
| rRNA Content                       | ?                    | Verify that excess ribosomal content is filtered/normalized across samples to ensure that alignment rates and subsequent normalization of data is not skewed.                                                                                                                                                                                                                                                |
| Transcript Coverage and 5’-3’ Bias | [Qualimap][qualimap] | Libraries prepared with polyA selection may have higher biased expression in 3’ region.  If reads primarily accumulate at the 3’ end of transcripts (in poly(A)-selected samples), this might indicate the starting RNA was of low quality.                                                                                                                                                                  |
| Junction Analysis                  | [Qualimap][qualimap] | Analysis of known, partly known, and novel junction positions in spliced alignments.                                                                                                                                                                                                                                                                                                                         |
| Strand Specificity                 | RSeQC                | Verification/sanity check of how reads were stranded for the RNA sequencing (stranded or unstranded protocol).                                                                                                                                                                                                                                                                                               |
| GC Content Bias                    | ?                    | GC profiles are typically remarkably stable.  Even small/minor deviations could indicate a problem with the library used (or bacterial contamination).                                                                                                                                                                                                                                                       |

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

These are generic instructions for running each of the tools in our pipeline.
We run our pipeline in a series of QC scripts that are tailored for our compute
cluster, so those commands may not apply elsewhere.  Instead we've supplied
examples of the commands used to each package.  Our default memory is 80G and we
employ 4 threads for these processes.

## Dependencies

We presume anaconda is available and installed. If not please follow the link to [anaconda](https://www.anaconda.com/) first.

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

## Workflow

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
| `qualimap rnaseq`        | Screen for RNA-seq bias and junction analysis              |
| `rseqc infer_experiment` | Determine RNA-SEQ strandedness and reads                   |
| `multiqc`                | Report aggregation                                         |

Note: Specific options such as memory size thresholds and thread count are below.

1. Run `samtools quickcheck` to ensure that input BAMs are relatively
   well-formed (for instance, to ensure a header and EOF marker exist).

   ```bash
   samtools quickcheck $BAM
   ```

2. Use Picard's `ValidateSamFile` tool to ensure the inner contents of the BAM
   file are well-formed.

   ```bash
   picard ValidateSamFile \
       I=$BAM \                                        # specify bam file
       MODE=SUMMARY\                                   # concise output
       INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE \   # lower stringency faster processing time
       IGNORE=INVALID_PLATFORM_VALUE                   # Validations to ignore.
   ```

3. Run `samtools flagstat` to gather general statistics such as alignment
   percentage.

   ```bash
   samtools flagstat $BAM 
   ```

4. Run `fastqc` to collect sequencing and library-related statistics. These are
   only for informational purposes — as stated above, we typically do not remove
   samples based on this information (with rare exception), as the
   sequencing-related QC work was done upstream in the genomics lab.

   ```bash
   fastqc $BAM
   ```

5. Run `qualimap bamqc` to gather more in-depth statistics about read stats,
   coverage, mapping quality, mismatches, etc.

   ```bash
   qualimap bamqc -bam $BAM \          # bam filename
       -nt $NUM_THREADS \              # threads requested
       -nw 400                         # number of windows
   ```

6. If RNA-seq data, run `qualimap rnaseq` to gather QC statistics that are
   tailored for RNA-seq files.

   ```bash
   qualimap rnaseq --java-mem-size=$MEM_SIZE \ # memory
       -bam $BAM \                             # bam filename
       -gtf $GTF_REF                           # transcript definition file
       -pe                                     # specify paired end if paired end
   ```

7. If RNA-seq, run `ngsderive strandedness` to determine a backwards-computed
   strandedness of the RNA-seq experiment.

   ```bash
   ngsderive strandedness
   ```

8. Compute the md5 checksum of the file.

   ```bash
   md5sum $BAM
   ```

9. Combine all of the above metrics using `multiqc`.

   ```bash
   multiqc . # recurse all files in '.'
   ```

# Items Still In-Progress

- [ ] Analysis tools for other types of sequencing (ChIP seq)
- [ ] Useful metadata from various stages (sample collection, laboratory, pre-sequencing, sequencing, post-sequencing)

# Outstanding Questions

- What thresholds or metrics differentiate a poor-quality sample from a high-quality one?
- What other metrics or properties would be valuable?
- What is best way to define and handle outliers?
- What is the best way to examine cohort integrity? This means experimental category-based tests of samples to find outliers that are of sufficient quality if examined alone.  Outliers in this case may indicate classification errors or rare biological conditions.  Which metrics are best tested here?

[samtools]: http://www.htslib.org/doc/samtools.html
[fastqc]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[qualimap]: http://qualimap.bioinfo.cipf.es/doc_html/command_line.html
[picard]: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_ValidateSamFile.php
[multiqc]: https://multiqc.info/