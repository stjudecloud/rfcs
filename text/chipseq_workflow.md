# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
- [Specification](#specification)

# Introduction

This RFC lays out the specification for the ChIP-Seq mapping pipeline.

# Motivation

To provide the epigenetics community access to data from ChIP-Seq experiments performed at St. Jude Children's Research Hospital, we propose the following data harmonization pipeline. The goal of this pipeline is to provide harmonized alignments for ChIP-Seq data. For this pipeline, we will make no recommendations on downstream analysis, focusing instead on harmonizing the underlying sequencing data and leaving analysis decisions to the user.

# Discussion

## Aligner choice

The epigenetics community and the results analyses of ChIP-seq experiments are highly varied. Unlike RNA-Seq or DNA-Seq, there is no standard mapping method for the analysis of human data. The community primarily uses either [Bowtie], [Bowtie 2], or [BWA] for alignment. These aligners use similar internal data structures for mapping via an FM-index.

To provide harmonized ChIP-Seq data in St. Jude Cloud, we need to select a single alignment method. Therefore we investigated how alignment was handled by other large projects.

[Bowtie]: http://bowtie-bio.sourceforge.net/index.shtml
[Bowtie 2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[BWA]: http://bio-bwa.sourceforge.net/

### Open source

Open source projects are ones that publish the source code of their analysis pipelines. This allows anyone to freely explore and review the details of each step.

| Project                                                   | Stable release     | Aligner                                                                                                                      |
| --------------------------------------------------------- | ------------------ | ---------------------------------------------------------------------------------------------------------------------------- |
| [ENCODE-DCC/chip-seq-pipeline] [[1][chip-seq-pipeline-1]] | 1.4 (2018-03-08)   | BWA 0.7.10 (2014-07-14) (bwa-aln)                                                                                            |
| [ENCODE-DCC/chip-seq-pipeline2]                           | 1.8.0 (2021-03-26) | Bowtie 2 2.3.4.3 (2018-09-17) (default)<br />BWA 0.7.17 (2017-11-07) (If read-length >= 70 bp, bwa-mem; otherwise, bwa-aln.) |
| [nf-core/chipseq] [[1][nf-core/chipseq-1]]                | 1.2.1 (2020-07-29) | BWA 0.7.17 (2017-11-07) (bwa-mem)                                                                                            |
| [stjude/SEASEQ]                                           | 1.1 (2021-04-19)   | Bowtie 1.2.3 (2019-07-05)                                                                                                    |

[ENCODE-DCC/chip-seq-pipeline]: https://github.com/ENCODE-DCC/chip-seq-pipeline
[chip-seq-pipeline-1]: http://doi.org/10.1101/gr.136184.111
[ENCODE-DCC/chip-seq-pipeline2]: https://github.com/ENCODE-DCC/chip-seq-pipeline2
[nf-core/chipseq]: https://nf-co.re/chipseq
[nf-core/chipseq-1]: https://www.nature.com/articles/s41587-020-0439-x
[stjude/SEASEQ]: https://github.com/stjude/seaseq

### Description only

Some projects only provide written descriptions of their analysis pipelines, either published in a paper or some other online document. These projects do not make their source code publicly available or easily accessible.

| Project                                             | Publication year | Aligner                     |
| --------------------------------------------------- | ---------------- | --------------------------- |
| [ChIP-Atlas] [[1][chip-atlas-1], [2][chip-atlas-2]] | 2018             | Bowtie 2 2.2.2 (2014-04-10) |
| [ChIPSummitDB] [[1][chipsummitdb-1]]                | 2020             | BWA                         |
| [Cistrome Project] [[1][cistrome-project-1]]        | 2018             | BWA                         |
| [Roadmap Epigenomics Project] [[1][rep-1]]          | 2015             | PASH 3.0 (2010)             |

[ChIP-Atlas]: https://chip-atlas.org/
[chip-atlas-1]: https://doi.org/10.15252/embr.201846255
[chip-atlas-2]: https://github.com/inutano/chip-atlas/wiki#2-primary-processing
[ChIPSummitDB]: http://summit.med.unideb.hu/summitdb/
[chipsummitdb-1]: https://doi.org/10.1093/database/baz141
[Cistrome Project]: http://cistrome.org/
[cistrome-project-1]: https://doi.org/10.1093/nar/gky1094
[Roadmap Epigenomics Project]: http://www.roadmapepigenomics.org/
[rep-1]: https://doi.org/10.1038/nature14248

Within St. Jude, Computational Biology has historically used BWA to map ChIP-Seq data. Also, the St. Jude Center for Applied Bioinformatics uses BWA as their standard aligner for ChIP-Seq experiments.

## Multiple mapped reads

To reduce false positives, ChIP-Seq experiments commonly filter out reads that are not uniquely mapped.  In order to simplify the ChIP-Seq workflow and to provide the most widely usable alignment data, however, we do not filter reads.

# Specification

## Required Metadata

Due to the nature of the experimental methods, ChIP-Seq requires additional metadata to be understood and used for analysis. The following information must be provided for all ChIP-Seq samples to enable their release on St. Jude Cloud.

* Antibody manufacturer/vendor
* Antibody product/catalog number
* Antibody lot number
* Antibody target/mark

## Dependencies

If you'd like the full `conda` environment, you can install it using the following command. Obviously, you'll need to install [anaconda](https://www.anaconda.com/) first.

```bash
conda create -n chipseq-mapping \
    -c conda-forge \
    -c bioconda \
    picard==2.20.2 \
    samtools==1.9 \
    bwa==0.7.17 \
    deeptools==3.5.0 \
    -y
```

Additionally, you will want to install our `fqlib` library to check that FastQ files are properly paired and have no duplicate read names. Installation of the [Rust](https://rustup.rs/) programming language is required.

```bash
cargo install --git https://github.com/stjude/fqlib.git --tag v0.3.1
```

## Reference files

The following reference files are used as the basis of the ChIP-Seq Workflow:

- Similarly to all analysis pipelines in St. Jude Cloud, we use the `GRCh38_no_alt` analysis set for our reference genome. You can get a copy of the file [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). Additionally, you can get the file by running the following commands:

  ```bash
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O GRCh38_no_alt.fa.gz
  gunzip GRCh38_no_alt.fa.gz

  echo "a6da8681616c05eb542f1d91606a7b2f  GRCh38_no_alt.fa" > GRCh38_no_alt.fa.md5
  md5sum -c GRCh38_no_alt.fa.md5
  # > GRCh38_no_alt.fa: OK
  ```

  - Last, the following command is used to prepare the BWA index file:

  ```bash
  bwa index GRCh38_no_alt.fa
  ```

## Workflow

Here are the resulting steps in the ChIP-Seq Workflow pipeline. There might be slight alterations in the actual implementation, which can be found in [the St. Jude Cloud workflows repository](https://github.com/stjudecloud/workflows/blob/master/workflows/chipseq/chipseq-standard.wdl).

1. Run `picard ValidateSam` on the incoming BAM to ensure that it is well-formed enough to convert back to FastQ.

   ```bash
   picard ValidateSamFile I=$INPUT_BAM \                # Input BAM.
                     IGNORE=INVALID_PLATFORM_VALUE \    # Validations to ignore.
                     IGNORE=MISSING_PLATFORM_VALUE
   ```

2. Run `samtools quickcheck` on incoming BAM to ensure it is properly formatted.

   ```bash
   samtools quickcheck $INPUT_BAM
   ```

3. Split BAM file into multiple BAMs on the different read groups using `samtools split`. See [the samtools documentation](http://www.htslib.org/doc/samtools.html) for more information.

   ```bash
   samtools split -u $UNACCOUNTED_BAM_NAME \ # Reads that do not belong to a read group or the read group is unrecognized go here.
                  -f '%*_%!.%.'              # Format of output BAM file names.
   ```

   If the BAM has unaccounted reads, those reads will need to be removed and the samtools split step will need to be rerun.

4. Run Picard `SamToFastq` on each of the BAMs generated in the previous step.

   ```bash
      picard SamToFastq \
             INPUT=$INPUT_BAM \
             FASTQ=$FASTQ_R1 \
             SECOND_END_FASTQ=$FASTQ_R2 \
             RE_REVERSE=true \
             VALIDATION_STRINGENCY=SILENT
   ```

5. Run `fq lint` on each of the FastQ pairs generated in the previous step as a quality check. You can see the checks that the `fq` tool performs [here](https://github.com/stjude/fqlib/blob/master/README.md#validators).

   ```bash
   fq lint $FASTQ_R1 $FASTQ_R2            # Files for read 1 and read 2. Read 2 is optional.
   ```

6. Run the `BWA` alignment algorithm.

   ```bash
   bwa aln -t ${ncpu} $INDEX_PREFIX \
        $ALL_FASTQ_R1 > sai_1 \           # FastQ files, separated by comma if there are multiple files.

   # If paired end data
   bwa aln -t ${ncpu} $INDEX_PREFIX \
        $ALL_FASTQ_R2 > sai_2 \           # FastQ files, separated by comma if there are multiple files.

   # For single end data
   bwa samse \
        -r $READ_GROUP_STRING \
        $INDEX_PREFIX \
        $SAI \
        $ALL_FASTQ_R1 | samtools view -hb --threads ${ncpu} > $BWA_BAM

   # For paired end data
   bwa sampe \
        -r $READ_GROUP_STRING \
        $INDEX_PREFIX \
        $SAI_R1 \
        $SAI_R2 \
        $ALL_FASTQ_R1 \
        $ALL_FASTQ_R2 | samtools view -hb --threads ${ncpu} > $BWA_BAM
   ```

7. Run `picard CleanSam` on the `BWA`-aligned BAM file. Fixes soft-clipping beyond the end-of-reference and sets MAPQ to 0 for unmapped reads.

   ```bash
   picard -Xmx${JAVA_HEAP_SIZE}g CleanSam \
      I=${BAM} \                                      # Input BAM.
      O=${OUTPUT_FILENAME}                            # Output cleaned BAM
   ```

8. Run `picard MergeSamFiles` on aligned, cleaned BAM files.

   ```bash
      picard -Xmx${JAVA_HEAP_SIZE}g MergeSamFiles \
         ${INPUT_ARG} \                               # Write an INPUT= argument for each input BAM file
         OUTPUT=${OUTPUT_NAME} \                      # Nmae for combined BAM file
         SORT_ORDER=${SORT_ORDER} \                   # Set sort order (default=coordinate)
         USE_THREADING=${THREADING} \                 # Use specified number of threads
         VALIDATION_STRINGENCY=SILENT                 # Ignore validation errors
   ```

9. Index the coordinate-sorted BAM file.

   ```bash
   samtools index $BWA_SORTED_BAM                     # BWA-aligned, coordinate-sorted BAM.
   ```

10. Run `picard ValidateSamFile` on the aligned and sorted BAM file.

   ```bash
   picard ValidateSamFile I=$BWA_SORTED_BAM \         # BWA-aligned, coordinate-sorted BAM.
                  IGNORE=INVALID_PLATFORM_VALUE \     # Validations to ignore.
                  IGNORE=MISSING_PLATFORM_VALUE
   ```

11. Run `bamCoverage` to generate bigwig file.

    ```bash
    bamCoverage --bam ${MERGED_BAM} \                 # Input BAM file
                --outFileName ${PREFIX}.bw \          # Output bigwig filename
                --outFileFormat bigwig \              # Set output format to bigwig
                --numberOfProcessors "max"            # Utilize all available processors
    ```
