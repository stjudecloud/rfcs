# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
- [Specification](#specification)
- [Appendix](#appendix)

# Introduction

This RFC lays out the specification for the ChIP-Seq mapping pipeline.

# Motivation

To provide the epigenetics community access to data from ChIP-Seq experiments performed at St. Jude Children's Research Hospital, we propose the following data harmonization pipeline. The goal of this pipeline is to provide harmonized alignments for ChIP-Seq data. For this pipeline, we will make no recommendations on downstream analysis, focusing instead on harmonizing the underlying sequencing data and leaving analysis decisions to the user.

# Discussion

## Aligner choice

The epigenetics community and the results analyses of ChIP-seq experiments are highly varied. Unlike RNA-Seq or DNA-Seq, there is no standard mapping method for the analysis of human data. The community primarily uses either `bowtie2` or `bwa` for alignment. Both of these aligners use similar internal data structures for mapping via an FM-index.

To provide harmonized ChIP-Seq data in St. Jude Cloud, we need to select a single alignment method. Therefore we investigated how alignment was handled by other large projects.

| Resource                    | Aligner  | Reference                                                                                             |
| --------------------------- | -------- | ----------------------------------------------------------------------------------------------------- |
| ENCODE                      | BWA      | https://genome.cshlp.org/content/22/9/1813.full, https://www.encodeproject.org/pipelines/ENCPL220NBH/ |
| ChIP-Atlas                  | Bowtie2  | https://chip-atlas.org/, https://github.com/inutano/chip-atlas/wiki#primary_processing_doc            |
| ChIPSummitDB                | BWA      | https://academic.oup.com/database/article/doi/10.1093/database/baz141/5700342#191914006               |
| Roadmap Epigenomics Project | Pash 3.0 | https://www.nature.com/articles/nature14248#Sec12                                                     |
| Cistrome                    | BWA      | http://cistrome.org/db/#/, https://academic.oup.com/nar/article/47/D1/D729/5193328#129642788          |

Within St. Jude, Computational Biology has historically used BWA to map ChIP-Seq data. Also, the St. Jude Center for Applied Bioinformatics uses BWA as their standard aligner for ChIP-Seq experiments.

## Multiple mapped reads

To reduce false positives, ChIP-Seq experiments commonly filter out reads that are not uniquely mapped.  In order to simplify the ChIP-Seq workflow and to provide the most widely usable alignment data, however, we do not filter reads.

# Specification

## Dependencies

If you'd like the full `conda` environment, you can install it using the following command. Obviously, you'll need to install [anaconda](https://www.anaconda.com/) first.

```bash
conda create -n chipseq-mapping \
    -c conda-forge \
    -c bioconda \
    picard==2.20.2 \
    samtools==1.9 \
    bwa==0.7.17 \
    ngsderive==1.0.2 \
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
2. Split BAM file into multiple BAMs on the different read groups using `samtools split`. See [the samtools documentation](http://www.htslib.org/doc/samtools.html) for more information.

   ```bash
   samtools split -u $UNACCOUNTED_BAM_NAME \ # Reads that do not belong to a read group or the read group is unrecognized go here.
                  -f '%*_%!.%.'              # Format of output BAM file names.
   ```

   If the BAM has unaccounted reads, those will need to be removed and the samtools split step will need to be rerun.

3. Run Picard `SamToFastq` on each of the BAMs generated in the previous step.

   ```bash
      picard SamToFastq \
             INPUT=$INPUT_BAM \
             FASTQ=$FASTQ_R1 \
             SECOND_END_FASTQ=$FASTQ_R2 \
             RE_REVERSE=true \
             VALIDATION_STRINGENCY=SILENT
   ```

4. Run `fq lint` on each of the FastQ pairs generated in the previous step as a quality check. You can see the checks that the `fq` tool performs [here](https://github.com/stjude/fqlib/blob/master/README.md#validators).

   ```bash
   fq lint $FASTQ_R1 $FASTQ_R2 # Files for read 1 and read 2.
   ```

5. Run the `BWA` alignment algorithm.

   ```bash
   bwa aln $INDEX_PREFIX \
        $ALL_FASTQ_R1 $ALL_FASTQ_READ2 \ # FastQ files, separated by comma if there are multiple files. The order of your R1 and R2 files must match!

   # For single end data
   bwa samse \
        $INDEX_PREFIX \
        $SAI \
        $ALL_FASTQ_R1 | samtools view -hb --threads ${ncpu} -o $BWA_BAM

   # For paired end data
   bwa sampe \
        $INDEX_PREFIX \
        $SAI_R1 \
        $SAI_R2 \
        $ALL_FASTQ_R1 \
        $ALL_FASTQ_READ2 | samtools view -hb --threads ${ncpu} -o $BWA_BAM
   ```

6. Run `picard SortSam` on the `BWA`-aligned BAM file.

   ```bash
   picard SortSam I=$BWA_BAM \                   # Input BAM.
                  O=$BWASORTED_BAM \             # Coordinate-sorted BAM.
                  SO="coordinate" \              # Specify the output should be coordinate-sorted
                  CREATE_INDEX=false \           # Explicitly do not create an index at this step, in case the default changes.
                  CREATE_MD5_FILE=false \        # Explicitly do not create an md5 checksum at this step, in case the default changes.
                  COMPRESSION_LEVEL=5 \          # Explicitly set the compression level to 5, although, at the time of writing, this is the default.
                  VALIDATION_STRINGENCY=SILENT   # Turn off validation stringency for this step.
   ```

7. Index the coordinate-sorted BAM file.

   ```bash
   samtools index $BWA_SORTED_BAM # BWA-aligned, coordinate-sorted BAM.
   ```

8. Run `picard ValidateSamFile` on the aligned and sorted BAM file.

   ```bash
   picard ValidateSamFile I=$BWA_SORTED_BAM \     # BWA-aligned, coordinate-sorted BAM.
                  IGNORE=INVALID_PLATFORM_VALUE \ # Validations to ignore.
                  IGNORE=MISSING_PLATFORM_VALUE
   ```

9. Run `picard MarkDuplicates` on the `BWA`-aligned BAM file.

   ```bash
   picard MarkDuplicates I=$BWA_SORTED_BAM \            # Input BAM.
                         O=$MARKED_BAM \                # Duplicate-marked output BAM.
                         VALIDATION_STRINGENCY=SILENT \ # Turn of validation stringency for this step.
                         CREATE_INDEX=false \           # Explicitly do not create an index at this step, in case the default changes.
                         CREATE_MD5_FILE=false \        # Explicitly do not create an md5 checksum at this step, in case the default changes.
                         COMPRESSION_LEVEL=5 \          # Explicitly set the compression level to 5, although, at the time of writing, this is the default.
                         METRICS_FILE=$METRICS_FILE \   # Location for the metrics file produced by MarkDuplicates.
   ```

10. Run `bamCoverage` to generate bigwig file.

    ```bash
    bamCoverage --bam ${MARKED_BAM} \              # Input BAM file
                --outFileName ${prefix}.bw \       # Output bigwig filename
                --outFileFormat bigwig \           # Set output format to bigwig
                --numberOfProcessors "max"         # Utilize all available processors
    ```

# Appendix
