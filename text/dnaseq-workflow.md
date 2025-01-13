# DNA-Seq Harmonization Pipeline

# Introduction

This RFC lays out the specification for a harmonization workflow for DNA-Seq (typically, Whole Genome (WGS) and Whole Exome (WES)) data. This will replace the current, blackbox offering from Microsoft Genomics that has now been deprecated.

# Motivation

## Microsoft Genomics
Previously, we have relied on the [Microsoft Genomics Service](https://www.microsoft.com/en-us/genomics/). Their implementation is a wrapper of the Burrows-Wheeler Aligner (BWA) and the Genome Analysis Toolkit (GATK). However, in late 2024, we were notified of Microsoft's intention to retire the product and stop supporting it. This, combined with difficulties getting support for the product and the obfuscation of the pipeline details, has created a need for our own implementation of a DNA-Seq pipeline.

## Genomes and gene annotations
In the past we have used multiple versions of the human genome. These have been primarily hg38 (GRCh38), but with different additional chromosomes (such as the Epstein-Barr Virus (EBV) contig). There have also been incremental improvements to the hg38 genome that have been released over the past years since its initial release. Additionally, since the initiation of this project, the Telomere-to-Telomere (T2T) consortium has announced and published the first gapless, complete human genome sequence.

# Discussion

# Specification

## Dependencies

## Reference files

## Workflow

1. Run `picard ValidateSam` on the incoming BAM to ensure that it is well-formed enough to convert bak to FastQ.

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

   If the BAM has unaccounted reads, those will need to be triaged and this step will need to be rerun.

3. Run Picard `SamToFastq` on each of the BAMs generated in the previous step.
   ```bash
      picard SamToFastq \
             INPUT=$INPUT_BAM \
             FASTQ=$FASTQ_R1 \
             SECOND_END_FASTQ=$FASTQ_R2 \
             RE_REVERSE=true \
             VALIDATION_STRINGENCY=SILENT
   ```

4. Run `fq lint` on each of the FastQ pairs that was generated in the previous step as a sanity check. You can see the checks that the `fq` tool performs [here](https://github.com/stjude/fqlib/blob/master/README.md#validators).

   ```bash
   fq lint $FASTQ_R1 $FASTQ_R2 # Files for read 1 and read 2.
   ```

5. Split FastQ pairs into smaller subsets (default: 10,000,000 reads) to enable additional parallelism and reduce turnaround time per-sample.

    ```bash
        let "lines = $READS_PER_FILE * 4"
        zcat ~{fastq} | split -l $lines -d -a 6 - $PREFIX

        for file in "$PREFIX"*; do
            mv "$file" "${file}.fastq"
            echo "gzip ${file}.fastq" >> cmds
        done

        parallel --jobs $NCPU < cmds
    ```

6. Run the `bwa` alignment algorithm.

    ```bash
        bwa mem \
            -t "$NCPU" \
            -R $READ_GROUP \
            bwa_db/"$PREFIX" \
            $READ_ONE \
            $READ_TWO \
            | samtools view --no-PG --threads "$NCPU" -hb - \
            > $OUTPUT_BAM
    ```

7. Sort individual BAM files with `picard`.

    ```bash
        picard -Xmx$JAVA_HEAP_SPACEg SortSam \
            -I $BAM \
            -O $OUTPUT_BAM \
            -SO "coordinate" \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --VALIDATION_STRINGENCY "SILENT"
    ```

8. Merge individual BAM files with `samtools`.

    ```bash
        samtools merge \
            --threads "$NCPU" \
            $PREFIX.bam \
            $BAMS
    ```

9.  Index final BAM file with `samtools`.

    ```bash
        samtools index --threads "$NCPU" $BAM $OUTPUT_BAM
    ```
