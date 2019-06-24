# Table of Contents <!-- omit in toc -->

- [Introduction](#Introduction)
- [Motivation](#Motivation)
- [Proposed changes](#Proposed-changes)
  - [Tool additions updates](#Tool-additions-updates)
  - [Update reference files.](#Update-reference-files)
  - [QC and quality of life improvements](#QC-and-quality-of-life-improvements)
  - [Various other changes](#Various-other-changes)
- [Summary](#Summary)
  - [Todo](#Todo)
  - [Outstanding Questions](#Outstanding-Questions)

# Introduction

This RFC lays out some thoughts I've been collecting about how to improve the RNA-Seq mapping pipeline based on (a) new version of tools/reference files and (b) feedback from the community.

# Motivation

* **Tool additions and updates.** The tools we use are woefully out of date (~2 years old). We should reap the benefits of new tools if possible. Additionally, there is some new functionality in the area of QC and validation that I'd like to add. See the [section below](#Tool-additions-updates) for more details.
  * Note that all of the tools used in the RNA-Seq Workflow v1.0 were the latest available version.
* **Updated reference files.** No changes have really been made to the `GRCh38_no_alt` analysis set FASTA. However, two major releases of the GENCODE gene model have transpired since v1.0 (we are now on [GENCODE v30](https://www.gencodegenes.org/human/release_30.html)).
* **QC and quality of life improvements based on feedback from the community.** Many interactions with the community have impacted the thoughts in this release. One of the most important themes in the RNA-Seq Workflow v2.0 proposal is the emphasis on QC and quality of life improvements (e.g. `fq lint`, generation and publication of md5sums).

# Proposed changes

## Tool additions updates

* Add `fq v0.2.0` ([Released](https://github.com/stjude/fqlib/releases/tag/v0.2.0) November 28, 2018) to be used in validation of the `picard SamToFastq` step.
* Add `qualimap v.2.2.2` ([Source](https://bitbucket.org/kokonech/qualimap/)) to use the `bamqc` and `rnaseq`tools. We've already been doing this in our QC pipeline, but I'd like to pull this out formally into the RNA-Seq workflow (mostly for documentation purposes).
* Update `STAR 2.5.3a` ([Released](https://github.com/alexdobin/STAR/releases/tag/2.5.3a) March 17, 2017) to `STAR 2.7.1a` ([Released](https://github.com/alexdobin/STAR/releases/tag/2.7.1a) May 15, 2019).
* Update `samtools 1.4.0` ([Released](https://github.com/samtools/samtools/releases/tag/1.4) March 13, 2017) to `samtools 1.9` ([Released](https://github.com/samtools/samtools/releases/tag/1.9) July 18, 2018).
  * Updating the samtools version whenever possible is of particular interest to me, due to the perceived fragility of the samtools code (although it has seemed to get better over the last year or so here).
* Update `picard 2.9.4` ([Released](https://github.com/broadinstitute/picard/releases/tag/2.9.4) June 15, 2017) to `picard 2.20.2` ([Released](https://github.com/broadinstitute/picard/releases/tag/2.20.2) May 28, 2019)

## Update reference files.

* `GENCODE v28` to `GENCODE v30`. Routine updates to the gene model from ENCODE.
* Previously, the we were filtering out anything not matching "level 1" or "level 2" from the gene model. This was due to best practices outlined during our RNA-Seq Workflow v1.0 discussions. I propose we revert this for the following reasons:
  * The first sentence in section 2.2.2 of the [STAR 2.7.1.a manual](https://github.com/alexdobin/STAR/blob/2.7.1a/doc/STARmanual.pdf): "The use of the most comprehensive annotations for a given species is strongly recommended". So it seems the author recommends you use the most comprehensive gene model.
  * Here is what [the GENCODE FAQ](https://www.gencodegenes.org/pages/faq.html) has to say about the level 3 annotations: "Ensembl loci where they are different from the Havana annotation or where no Havana annotation can be found". Given that the GENCODE geneset is the union of automated annotations from the `Ensembl-genebuild` and manual curation of the `Ensembl-Havana` team, this level should be harmless in the event that levels 1 & 2 don't apply.
  * Last, the [GDC documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#rna-seq-alignment-command-line-parameters) does not currently describe any filtering of their GENCODE v22 GTF (see the section labeled "Step 1: Building the STAR index").

## QC and quality of life improvements

* Add `picard ValidateSamFile` to the checks after the `STAR` alignment and `picard MarkDuplicates` steps. The criticism internally is that `ValidateSamFile` is quite stringent and often errors with concerns we don't care about. I'm testing this out as I develop the pipeline, and so far, I've found the following warnings to be ignore-worthy:
  * `INVALID_PLATFORM_VALUE` is pretty annoying. It just complains if a read group doesn't contain a `PL` attribute. I'm not sure it's worth going back and fixing these.
* For dependency management, I'd like to propose we move to using `conda` until we are decided on which workflow language we will support and until we get further down the road of building standard docker images. All packages should be available within the `defaults`, `conda-forge`, and `bioconda` repositories.
* Add a checksum algorithm and publish the results in the data browser. Currently, I'm proposing we generate the `md5sum` checksum. However, we should consider the use of a non-broken hashing algorithm (see [the related question below](#Outstanding-Questions)).

## Various other changes

* Removed a section of the pipeline that reformatted the header of the BAM to be cleaner. `STAR` outputs a header that is formatted fine already, and I found this code to just be an area where an error could be introduced for little benefit.
* I removed a section of custom code that checks for duplicate read groups. `picard ValidateSamFile` does this for you (see [the documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=7571) for this tool. Specifically, the `DUPLICATE_READ_GROUP_ID` error).


# Summary

The following reference files are used as the basis of the RNA-Seq Workflow v2.0:

* Similarly to all analysis pipelines in St. Jude Cloud, we use the `GRCh38_no_alt` analysis set for our reference genome. You can get a copy of the file [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). Additionally, you can get the file by running the following commands:

   ```bash
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O GRCh38_no_alt.fa.gz
   gunzip GRCh38_no_alt.fa.gz

   echo "a6da8681616c05eb542f1d91606a7b2f  GRCh38_no_alt.fa" > GRCh38_no_alt.fa.md5
   md5sum -c GRCh38_no_alt.fa.md5
   # > GRCh38_no_alt.fa: OK
   ```

* For the gene model, we use the GENCODE v30 "comprehensive gene annotation" GTF for the "CHR" regions. You can get a copy of the gene annotation file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz). Additionally, you can get the file by running the following commands:

   ```bash
   wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz
   gunzip gencode.v30.annotation.gtf.gz

   echo "63770a3d2c6adb4d9d1bbc9ba3bd4adf  gencode.v30.annotation.gtf" > gencode.v30.annotation.gtf.md5
   md5sum -c gencode.v30.annotation.gtf.md5
   # > gencode.v30.annotation.gtf: OK
   ```

If you'd like the full `conda` environment, you can install it using the following command. Obviously, you'll need install [anaconda](https://www.anaconda.com/) first.

```bash
conda create -n star-mapping \
    -c conda-forge \
    -c bioconda \
    picard==2.20.2 \
    samtools==1.9 \
    star==2.7.1a \
    qualimap==2.2.2c \
    multiqc==1.7 \
    -y
```


Here are the resulting steps in the RNA-Seq Workflow v2.0 pipeline.

1. Prepare the STAR index file.

   ```bash
   STAR --runMode genomeGenerate \                    # Use genome generation runMode.
        --genomeDir $OUTPUT_DIR \                     # Specify an output directory.
        --runThreadN $NCPU \                          # Number of threads to use to build genome database.
        --genomeFastaFiles $FASTA \                   # A path to the GRCh38_no_alt.fa FASTA file.
        --sjdbGTFfile $GENCODE_GTF_V30 \              # GENCODE v30 gene model file, unmodified.
        --sjdbOverhang 125                            # Splice junction database overhang parameter, the optimal value is (Max length of RNA-Seq read-1).
   ```

2. Run `samtools quickcheck` on the incoming BAM to ensure that it is well-formed enough to convert back to FastQ.
3. Split BAM file into multiple BAMs on the different read groups using `samtools split`. See [the samtools documentation](http://www.htslib.org/doc/samtools.html) for more information.

    ```bash
    samtools split -u $UNACCOUNTED_BAM_NAME \ # Reads that do not belong to a read group or the read group is unrecognized go here.
                   -f '%*_%!.%.'              # Format of output BAM file names.
    ```

    If the BAM has unaccounted reads, those will need to be triaged and this step will need to be rerun.
4. Run `fq lint` on each of the FastQ pairs that was generated from step 3 as a sanity check. You can see the checks that the `fq` tool performs [here](https://github.com/stjude/fqlib/blob/master/README.md#validators)

    ```bash
    fq lint $FASTQ_R1 $FASTQ_R2 # Files for read 1 and read 2.
    ```
5. Run the `STAR` alignment algorithm.

    ```bash
    STAR --readFilesIn $ALL_READ1 $ALL_READ2 \     # FastQ files, separated by comma if there are multiple. The order of your R1 and R2 files has to match!
         --outSAMattrRGline $ALL_RG \              # Read group lines in the same order as `readFilesIn` (derived from earlier `samtools split` step).
         --genomeDir $STARDB \                     # Directory containing the STAR genome
         --runThreadN $NCPU \                      # Number of threads to use. You must request the correct amount from HPCF first!
         --outSAMunmapped Within \                 # Keep unmapped reads in the final BAM file.
         --outSAMstrandField intronMotif \         # Preserve compatibility with Cufflinks by including the XS attribute (strand derived from intron motif).
         --outSAMtype BAM SortedByCoordinate \     # Output a BAM file that is coordinate sorted.
         --outSAMattributes NH HI AS nM NM MD XS \ # Recommended SAM attributes to include for compatibility. Refer to manual for specifics.
         --outFilterMultimapScoreRange 1 \         # Ensures that all multi-mapped reads will need to share the mapping score.
         --outFilterMultimapNmax 20 \              # Max number of multi-mapped SAM entries for each read.
         --outFilterMismatchNmax 10 \              # Max number of mismatches allowed from alignment.
         --alignIntronMax 500000 \                 # Max intron size considered.
         --alignMatesGapMax 1000000 \              # Max gap allowed between two mates.
         --sjdbScore 2 \                           # Additional weight given to alignments that cross database junctions when scoring alignments.
         --alignSJDBoverhangMin 1 \                # Minimum overhang required on each side of a splice junction. Here, we require only one base on each side.
         --outFilterMatchNminOverLread 0.66 \      # 66% of the read must be perfectly matched to the reference sequence.
         --outFilterScoreMinOverLread 0.66 \       # Score must be greater than 66% of the read length. So for RL=100, the alignment must have a score > 66.
         --limitBAMsortRAM $RAM_LIMIT \            # Amount of RAM to use for sorting. Recommended value is [Max amount of RAM] - 5GB.
         --outFileNamePrefix $OUT_FILE_PREFIX      # All output files will have this path prepended.
         --twopassMode basic                       # Use STAR two-pass mapping technique (refer to manual).
    ```

6. Run `picard MarkDuplicates` on the `STAR`-aligned BAM file.
   
   ```bash
   picard MarkDuplicates I=$STAR_BAM \                  # Input BAM.
                         O=$MARKED_BAM \                # Duplicate-marked output BAM.
                         VALIDATION_STRINGENCY=SILENT \ # Turn of validation stringency for this step.
                         CREATE_INDEX=false \           # Explicitly do not create an index at this step, in case the default changes.
                         CREATE_MD5_FILE=false \        # Explicity do not create an md5 checksum at this step, in case the default changes.
                         COMPRESSION_LEVEL=5 \          # Explicitly set the compression level to 5, although, at the time of writing, this is the default.
                         METRICS_FILE=$METRICS_FILE \   # Location for the metrics file produced by MarkDuplicates.
   ```
7. Run `picard ValidateSamFile` on the aligned and marked BAM file.
   
   ```bash
   picard ValidateSamFile I=$INPUT_BAM \                # Input BAM.
                          IGNORE=INVALID_PLATFORM_VALUE # Validations to ignore.
   ```

8. Run `qualimap bamqc` and `qualimap rnaseq` QC for assistance in post-processing QC.

    ```bash
    qualimap bamqc -bam $INPUT_BAM \     # Input BAM. 
                   -outdir $OUTPUT_DIR \ # Output directory.
                   -nt $NCPUS \          # Number of CPUs to use.
                   -pe                   # All RNA-Seq data in St. Jude Cloud is currently paired end.
    ```

    and

    ```bash
    qualimap rnaseq -bam $INPUT_BAM \       # Input BAM.
                    -gtf $GENCODE_GTF_V30 \ # GENCODE v30 gene model file, unmodified.
                    -outdir $OUTPUT_DIR \   # Output directory.
                    -pe                     # All RNA-Seq data in St. Jude Cloud is currently paired-end.
    ```
9. Finally, generate the remaining files generally desired as output for the RNA-Seq Workflow.

   ```bash
   samtools flagstat $INPUT_BAM
   samtools index $INPUT_BAM
   md5sum $INPUT_BAM
   ```
10. Run `multiqc` across the following files for all samples in the cohort:

    * `fastqc`
    * `STAR`
    * `picard MarkDuplicates` and `picard ValidateSamFile`
    * `qualimap bamqc` and `qualimap rnaseq`
    * `samtools flagstat`

## Todo

- [ ] Investigation of impact for using ENCODE annotations post `GRCh38.p0` with the no alt analysis set. To measure this, I will see how many genes in the GENCODE gene model overlap with regions that are impacted by patches to the `GRCh38` genome.
- [ ] Is it a good idea/good investment of effort to remove absolute paths from the headers and leave just a relative path behind? For example, I think it would clean the header up significantly to change `/research/rgs01/project_space/zhanggrp/SJCloud/common/DataPreparation/RNA-Seq/PCGP-more-memory/data/SJAMLM7060_D/Output/SJAMLM7060_D.bam` to `/XXX/data/SJAMLM7060_D/Output/SJAMLM7060_D.bam`.
- [ ] Any read groups with `N/A` in the read group ID will cause `samtools split` to error out and try to create a file within a subdirectory. I'm considering functionality that will automatically replace any `N/A` string in a read group tag to `NA`.
- [x] Add `multiqc` to aggregate QC results.
- [x] Pin `qualimap` version.
- [ ] Update internal "STAR Best Practices" documentation.
- [ ] Update internal "Genome Data Files and Configuration" documentation.
- [ ] Update external documentation for RNA-Seq pipeline. Potentially break out the DNA-Seq and RNA-Seq workflows into their own file.
- [ ] Index files internally in TARTAn for `GRCh38_no_alt`.


## Outstanding Questions

* Any parameters we want to change during the STAR alignment step? I don't expect any, but we should explicitly discuss.
* Any additional end-to-end tests we'd like to add to the pipeline? Some have been suggested in the area of making sure certain transcripts quantify out the other side for a particular diagnosis.
* Should we include count quantification using `htseq-count` in the pipeline?
* This has just been an outstanding question of mine for a while — how big of an impact (if any at all) does the mismatch between the patch builds have? GENCODE v30 is built against `GRCh38.p12`, but obviously the no alt analysis set is derived from `GRCh38.p0` (with the PAR hard masked, the EBV chromosome added, etc.)?
* Should we be using `sha256` instead of `md5`? Just seems like using a non-broken hash algorithm would make sense. However, I'm not sure whether 
  * the `sha256sum` tool is sufficiently widespread enough, and
  * the benefit is worth the cost of breaking from the current community norm. However, this could also be a good thing for us to be forward thinking.