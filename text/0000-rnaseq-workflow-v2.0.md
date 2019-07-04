# Table of Contents <!-- omit in toc -->

- [Introduction](#Introduction)
- [Motivation](#Motivation)
- [Propsed Changes and Discussion](#Propsed-Changes-and-Discussion)
  - [Genomic tools](#Genomic-tools)
  - [Reference files](#Reference-files)
    - [GENCODE compatability](#GENCODE-compatability)
    - [Continued exclusion of automatically annotated features](#Continued-exclusion-of-automatically-annotated-features)
  - [QC and quality of life improvements](#QC-and-quality-of-life-improvements)
  - [Various other changes](#Various-other-changes)
- [Workflow description](#Workflow-description)
  - [Dependencies](#Dependencies)
  - [Reference files](#Reference-files-1)
  - [Workflow](#Workflow)
- [Items still in-progress](#Items-still-in-progress)
- [Outstanding questions](#Outstanding-questions)
- [Appendix](#Appendix)

# Introduction

This RFC lays out some thoughts I've been collecting about how to improve the RNA-Seq mapping pipeline based on (a) new version of tools/reference files and (b) feedback from the community.

# Motivation

* **Tool additions and updates.** The tools we use are woefully out of date (2 years old). We should reap the benefits of new tools if possible. Additionally, there is some new functionality in the area of QC and validation that I'd like to add. See the [section below](#Tool-additions-updates) for more details.
  * Note that all of the tools used in the RNA-Seq Workflow v1.0 were the latest available version.
* **Updated reference files.** No changes have really been made to the `GRCh38_no_alt` analysis set FASTA. However, three major releases of the GENCODE gene model have transpired since we released the first revision of the RNA-Seq workflow ([GENCODE v31](https://www.gencodegenes.org/human/release_31.html) is now out).
* **QC and quality of life improvements based on feedback from the community.** Many interactions with the community have impacted the thoughts in this release: 
  * A primary driver for the rewrite of the pipeline is the feedback we heard about the `ERCC SpikeIn` sequences. 
    * Popular tools such as `GATK` and `picard` are generally unhappy if the sequence dictionaries don't match perfectly. 
    * The inclusion of the `ERCC` genome in all of our RNA-Seq samples was violating that practice and causing problems for cross-target analyses.
    * Last, many of our samples do not currently contain the `ERCC` sequences.
    * After some discussion internally, we decided the best thing to do was to remove the ERCC genome by default. We are considering providing an ERCC version of the BAM for samples containing these sequences, but there is no consensus on whether it's worth it yet.
  * One of the most important themes in the RNA-Seq Workflow v2.0 proposal is the emphasis on QC and quality of life improvements (e.g. `fq lint`, generation and publication of md5sums).

# Propsed Changes and Discussion

## Genomic tools

* `fq v0.2.0` ([Released](https://github.com/stjude/fqlib/releases/tag/v0.2.0) November 28, 2018).
  * This tool will be used to validate the output of `picard SamToFastq`. `picard SamToFastq` does not currently catch all of the errors we wish to catch at this stage (such as duplicate read names in the FastQ file). Thus, we will leverage this tool to independently validate that the data is well-formed by our definition of that phrase.
* `rseqc v3.0.0` ([Source](http://rseqc.sourceforge.net/#download-rseqc)) will be added.
  * We have started using `infer_experiment.py` to infer strandedness from the data and ensure that the data matches what information we get from the lab.
* Added `qualimap v.2.2.2` ([Source](https://bitbucket.org/kokonech/qualimap/)).
  * Although we have been using `qualimap` quite heavily in our QC pipeline, we are formally adding this to the end of the RNA-Seq alignment workflow. The `bamqc` and `rnaseq` subcommands are both used.
* Update `STAR 2.5.3a` ([Released](https://github.com/alexdobin/STAR/releases/tag/2.5.3a) March 17, 2017) to `STAR 2.7.1a` ([Released](https://github.com/alexdobin/STAR/releases/tag/2.7.1a) May 15, 2019).
  * Upgraded to receive the benefits of bug fixes and software optimizations.
* Update `samtools 1.4.0` ([Released](https://github.com/samtools/samtools/releases/tag/1.4) March 13, 2017) to `samtools 1.9` ([Released](https://github.com/samtools/samtools/releases/tag/1.9) July 18, 2018).
  * Updating the samtools version whenever possible is of particular interest to me due to my perceived fragility of the samtools code (although it has seemed to get better over the last year or so here).
* Update `picard 2.9.4` ([Released](https://github.com/broadinstitute/picard/releases/tag/2.9.4) June 15, 2017) to `picard 2.20.2` ([Released](https://github.com/broadinstitute/picard/releases/tag/2.20.2) May 28, 2019).
  * Upgraded to receive the benefits of bug fixes and software optimizations.

## Reference files

### GENCODE compatability

One of the major discussions for this new version of the RNA-Seq workflow was the compatability of the `GRCh38_no_alt` reference genome with the latest `GENCODE` gene set. It was posed as the following question when the RFC was started:

    This has just been an outstanding question of mine for a while — how big of an impact (if any at all) does the mismatch between the patch builds have? GENCODE v31 is built against `GRCh38.p12` (see the in the title on the [webpage](https://www.gencodegenes.org/human/release_31.html)), but obviously the no alt analysis set is derived from `GRCh38.p0` (with the pseudoautosomal regions hard masked, the EBV chromosome added, etc.)
    
As is apparent in the question, the `GRCh38` reference genome is regularly patched with backwards compatable (with respect to coordinates) changes to the genome. Thus, gene models based on newer patches of `GRCh38` will not perfectly match the underlying nucleotide sequence in our reference genome wherever patches were applied.

As a first step, we researched what the best practices in the community are:

* The GDC are using GENCODE v22 against the GRCh38 no alt analysis set (sources for [gene model](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#rna-seq-alignment-command-line-parameters) and [reference genome](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files)). Although only separated by one patch, these two references appear to be discordant.
* ENCODE is using GENCODE v24 against the GRCh38 no alt analysis set ([source](https://www.encodeproject.org/data-standards/reference-sequences/)).

I reached out to the author of STAR, Alex Dobin, to get his opinion. You can read my question and his reply [here](https://github.com/alexdobin/STAR/issues/673).

Based on internal discussions and Dr. Dobin's reply, we considered three possible options:

1. We could try using the reference FASTA supplied with the respective GENCODE release as suggested by Dr. Dobin.
   * This was the most undesirable approach for a few reasons. The most prominent reason was that this reference genome did not immediately mirror the sequence dictionary or characteristics of the reference genome we use in our DNA-Seq pipeline, which was a major motiviation for this RFC. It would require a large amount of postprocessing of the GENCODE reference genome to convert to apply all of the no alt analysis set changes (e.g. converting to UCSC names, masking regions of the genome, inserting the EBV chromosome). This could leave room for the introduction of lots of strange errors, and there was no interest in getting into the business of genome generation (there is a reason they don't apply patches to the no alt analysis set).
2. We could downgrade the GENCODE gene model to be consistent with the no alt analysis patch ([v21](https://www.gencodegenes.org/human/release_21.html] would the correct version to use in this case).
   * The concordance of the reference sequences obviously made this choice an attractive option. However, `GENCODE v21` was released over 5 years ago (06.2014) and there were many valuable updates over that time. In particular, a quick look showed that there were many more transcripts and that the number of lncRNAs more than doubled (see appendix). We did not want to lose all of this forward progress if we could avoid it.
3. We could use the latest GENCODE release despite the mismatches between patches.
   * Given that there was no perfect option, this is the route we decided to take. The general consensus was that, after filtering the gene model for only known features (`level 1` and `level 2`), the effect of this discordance would be small enough to tolerate so that we could gain all of the knowledge accumulated since the older release. To quantify this, we measured the differences in gene expression (as measured by the `R^2` value between `GENCODE v21` and `GENCODE v31`) and splice junction detection (as measured by the number of splice junctions detected and the relative proportions of novel/partially novel/known junctions). In the current version of this RFC, I have not outlined the results, but I plan to in the future (see TODOs).

### Continued exclusion of automatically annotated features

Originally, I had posed this question to the group:

  * Previously, we were filtering out anything not matching "level 1" or "level 2" from the gene model. This was due to best practices outlined during our RNA-Seq Workflow v1.0 discussions. I propose we revert this for the following reasons:
    * The first sentence in section 2.2.2 of the [STAR 2.7.1.a manual](https://github.com/alexdobin/STAR/blob/2.7.1a/doc/STARmanual.pdf): "The use of the most comprehensive annotations for a given species is strongly recommended". So it seems the author recommends you use the most comprehensive gene model.
    * Here is what [the GENCODE FAQ](https://www.gencodegenes.org/pages/faq.html) has to say about the level 3 annotations: "Ensembl loci where they are different from the Havana annotation or where no Havana annotation can be found". Given that the GENCODE geneset is the union of automated annotations from the `Ensembl-genebuild` and manual curation of the `Ensembl-Havana` team, this level should be harmless in the event that levels 1 & 2 don't apply.
    * Last, the [GDC documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#rna-seq-alignment-command-line-parameters) does not currently describe any filtering of their GENCODE v22 GTF (see the section labeled "Step 1: Building the STAR index").

After discussion internally, we decided to continue removing `level 3` annotations by default. `level 3` features were considered too experimental to apply in the majority of use cases our end users explore. The small number of situations were one would want to include these in the counts weighed against the risk of false-positive discoveries in unverified genes + miscounted evidence that could have been counted towards verified features made it undesirable to include. 
   * For any users interested in these `level 3` features, you should rerun the alignment and quantification with these features included in your GTF. 

## QC and quality of life improvements

* Add `picard ValidateSamFile` to the checks after the `STAR` alignment and `picard MarkDuplicates` steps. The criticism internally is that `ValidateSamFile` is quite stringent and often errors with concerns we don't care about. I'm testing this out as I develop the pipeline, and so far, I've found the following warnings to be ignore-worthy:
  * `INVALID_PLATFORM_VALUE` is pretty annoying. It just complains if a read group doesn't contain a `PL` attribute. I'm not sure it's worth going back and fixing these.
* For dependency management, I'd like to propose we move to using `conda` until we are decided on which workflow language we will support and until we get further down the road of building standard docker images. All packages should be available within the `defaults`, `conda-forge`, and `bioconda` repositories.
* Add a checksum algorithm and publish the results in the data browser. Currently, I'm proposing we generate the `md5sum` checksum. However, we should consider the use of a non-broken hashing algorithm (see [the related question below](#Outstanding-Questions)).

## Various other changes

* Removed a section of the pipeline that reformatted the header of the BAM to be cleaner. `STAR` outputs a header that is formatted fine already, and I found this code to just be an area where an error could be introduced for little benefit.
* I removed a section of custom code that checks for duplicate read groups. `picard ValidateSamFile` does this for you (see [the documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=7571) for this tool. Specifically, the `DUPLICATE_READ_GROUP_ID` error).

# Workflow description

## Dependencies 

If you'd like the full `conda` environment, you can install it using the following command. Obviously, you'll need to install [anaconda](https://www.anaconda.com/) first.

```bash
conda create -n star-mapping \
    -c conda-forge \
    -c bioconda \
    picard==2.20.2 \
    samtools==1.9 \
    star==2.7.1a \
    qualimap==2.2.2c \
    multiqc==1.7 \
    rseqc==3.0.0 \
    fastqc==0.11.8-1 \
    -y
```

Additionally, you will want to install our `fqlib` library to check that FastQ files are properly paired and have no duplicate read names. Installation of the [Rust](https://rustup.rs/) programming language is required.

```bash
cargo install --git https://github.com/stjude/fqlib.git
```

## Reference files

The following reference files are used as the basis of the RNA-Seq Workflow v2.0:

* Similarly to all analysis pipelines in St. Jude Cloud, we use the `GRCh38_no_alt` analysis set for our reference genome. You can get a copy of the file [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). Additionally, you can get the file by running the following commands:

   ```bash
   wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O GRCh38_no_alt.fa.gz
   gunzip GRCh38_no_alt.fa.gz

   echo "a6da8681616c05eb542f1d91606a7b2f  GRCh38_no_alt.fa" > GRCh38_no_alt.fa.md5
   md5sum -c GRCh38_no_alt.fa.md5
   # > GRCh38_no_alt.fa: OK
   ```

* For the gene model, we use the GENCODE v31 "comprehensive gene annotation" GTF for the "CHR" regions and modify the file to include only "level 1" and "level 2" features (see discussion on why in the rationale below). You can get a copy of the gene annotation file [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz). For the exact steps to generate the gene model we use, you can run the following commands:

   ```bash
   wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
   echo "0aac86e8a6c9e5fa5afe2a286325da81  gencode.v31.annotation.gtf.gz" > gencode.v31.annotation.gtf.gz.md5
   md5sum -c gencode.v31.annotation.gtf.gz.md5
   # > gencode.v31.annotation.gtf.gz: OK

   gunzip gencode.v31.annotation.gtf.gz
   echo "4e22351ae216e72aa57cd6d6011960f8  gencode.v31.annotation.gtf" > gencode.v31.annotation.gtf.md5
   md5sum -c gencode.v31.annotation.gtf.md5
   # > gencode.v31.annotation.gtf: OK

   awk '/level 1/ || /level 2/' gencode.v31.annotation.gtf > gencode.v31.annotation.knownloci.gtf
   echo "8e5745b5b9b372936fb68858e89a0495  gencode.v31.annotation.knownloci.gtf" > gencode.v31.annotation.knownloci.gtf.md5
   md5sum -c gencode.v31.annotation.knownloci.gtf.md5
   # > gencode.v31.annotation.knownloci.gtf: OK
   ```

## Workflow

Here are the resulting steps in the RNA-Seq Workflow v2.0 pipeline.

1. Prepare the STAR index file.

   ```bash
   STAR --runMode genomeGenerate \                    # Use genome generation runMode.
        --genomeDir $OUTPUT_DIR \                     # Specify an output directory.
        --runThreadN $NCPU \                          # Number of threads to use to build genome database.
        --genomeFastaFiles $FASTA \                   # A path to the GRCh38_no_alt.fa FASTA file.
        --sjdbGTFfile $GENCODE_KNOWNLOCI_GTF_31 \     # GENCODE v31 gene model file including only level 1 and level 2 features (see generation steps above).
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
         --outFileNamePrefix $OUT_FILE_PREFIX \    # All output files will have this path prepended.
         --twopassMode Basic                       # Use STAR two-pass mapping technique (refer to manual).
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
8. Run `rseqc`'s `infer_experiment.py` to confirm that the lab's information on strandedness reflects what was is computed. Manually triage any descrepencies. This is particularly useful for historical samples.
    ```bash
    infer_experiment.py -i $INPUT_BAM -r $GENCODE_KNOWNLOCI_GTF_V31
    
    # Custom script to triage the following (these might be able to be simplified or improved, it's just my first stab):
    #   - If proportion of forward orientation evidence fraction is >= 0.8, assign "strand-specific-forward".
    #   - If proportion of reverse orientation evidence fraction is >= 0.8, assign "strand-specific-reverse".
    #   - If both proportions are between 0.6 and 0.4, assign "non-strand-specific".
    #   - Else flag for manual triage.
    ```

9. Run `qualimap bamqc` and `qualimap rnaseq` QC for assistance in post-processing QC. Note that for the `rnaseq` tool, we will need to include the strandedness for best results. The value received from the lab can generally be confirmed by the `infer_experiment.py` step above.

    ```bash
    qualimap bamqc -bam $INPUT_BAM \     # Input BAM. 
                   -outdir $OUTPUT_DIR \ # Output directory.
                   -nt $NCPUS \          # Number of CPUs to use.
                   -pe                   # All RNA-Seq data in St. Jude Cloud is currently paired end.
    ```

    and

    ```bash
    qualimap rnaseq -bam $INPUT_BAM \                  # Input BAM.
                    -gtf $GENCODE_KNOWNLOCI_GTF_V31 \  # GENCODE v31 gene model file including only level 1 and level 2 features (see generation steps above).
                    -outdir $OUTPUT_DIR \              # Output directory.
                    -oc qualimap_counts.txt \          # Counts as calculated by qualimap.
                    -p $COMPUTED \                     # Strandedness as specified by the lab and confirmed by "infer_experiment.py" above. Typically "strand-specific-reverse" for St. Jude Cloud data.
                    -pe                                # All RNA-Seq data in St. Jude Cloud is currently paired-end.
    ```
10. Next, `htseq-count` is run for the final counts file to be delivered:
    ```bash
    htseq-count -f bam \                   # Specify input file as BAM.
      -r pos \                             # Specify the BAM is position sorted.
      -s $COMPUTED \                       # Strandedness as specified by the lab and confirmed by "infer_experiment.py" above. Typically "reverse" for St. Jude Cloud data.
      -m union \                           # For reference, GDC uses "intersection-nonempty". Needs input from reviewers.
      -i gene_name \                       # I'd like to use the colloquial gene name here. For reference, GDC uses "gene_id" here. Needs input from reviewers.
      --secondary-alignments ignore \      # Elect to ignore secondary alignments. Needs input from reviewers.
      --supplementary-alignments ignore \  # Elect to ignore supplementary alignments. Needs input from reviewers.
      $INPUT_BAM                           # Input BAM file.
    ```
11. Generate the remaining files generally desired as output for the RNA-Seq Workflow.
    ```bash
    samtools flagstat $INPUT_BAM
    samtools index $INPUT_BAM
    md5sum $INPUT_BAM
    fastqc -f bam \     # Specify that we are working on a BAM file.
           -o $OUTDIR \ # Specify an out directory.
           -t $NCPU \   # Specify number of threads.
           $INPUT_BAM   # Input BAM we are QC'ing.
    ```
12. Run `multiqc` across the following files for all samples in the cohort:

    * `fastqc`
    * `STAR`
    * `picard MarkDuplicates` and `picard ValidateSamFile`
    * `qualimap bamqc` and `qualimap rnaseq`
    * `samtools flagstat`

# Items still in-progress

- [ ] Any read groups with `N/A` in the read group ID will cause `samtools split` to error out and try to create a file within a subdirectory. I'm considering functionality that will automatically replace any `N/A` string in a read group tag to `NA`.
- [x] Add `multiqc` to aggregate QC results.
- [x] Pin `qualimap` version.
- [x] Pin `fastqc` and add steps.
- [x] Pin `rseqc` and add steps.
- [ ] Pin `htseq-count` and add steps for read quantification.
- [ ] Update internal "STAR Best Practices" documentation.
- [ ] Update internal "Genome Data Files and Configuration" documentation.
- [ ] Update external documentation for RNA-Seq pipeline. Potentially break out the DNA-Seq and RNA-Seq workflows into their own file.
- [ ] Index files internally in TARTAn for `GRCh38_no_alt`.
- [ ] Add details about analysis done to choose v31 of the ENCODE gene model over v21.

# Outstanding questions

* Any parameters we want to change during the STAR alignment step? I don't expect any, but we should explicitly discuss.
* There are several parameters during the `htseq-count` step that need input from experts on what a sane default would be for our pipeline.
* I'd like to incorporate some manual end-to-end tests for our pipeline to evaluate results when we change parameters? Any ideas on this? The areas I have been thinking about are (a) gene expression quantification, (b) gene fusion detection, and (c) novel splice junction/isoform detection.
* Should we consider adding `qualimap counts` to QC counts across a cohort?
* Should we consider the generation and vending of `.bw` files by default?
* Should we be using `sha256` instead of `md5`? Just seems like using a non-broken hash algorithm would make sense. However, I'm not sure whether 
  * the `sha256sum` tool is sufficiently widespread enough, and
  * the benefit is worth the cost of breaking from the current community norm. However, this could also be a good thing for us to be forward thinking.

# Appendix

* Commands that were used in the comparison of feature types between `GENCODE v21` and `GENCODE v31`:

    ```bash
    gawk '$0 !~ /^#/ { features[$3] += 1; total += 1 } END {for (key in features) { print "Feature " key ": " features[key] " (" (features[key]/total)*100 "%)" }; print "Total: " total}' gencode.v21.annotation.gtf 
    # Feature exon: 1162114 (45.6341%)
    # Feature CDS: 696929 (27.3671%)
    # Feature UTR: 275100 (10.8027%)
    # Feature gene: 60155 (2.36217%)
    # Feature start_codon: 81862 (3.21457%)
    # Feature Selenocysteine: 114 (0.00447657%)
    # Feature stop_codon: 73993 (2.90557%)
    # Feature transcript: 196327 (7.7094%)
    # Total: 25465949
  
    gawk '$0 !~ /^#/ { features[$3] += 1; total += 1 } END {for (key in features) { print "Feature " key ": " features[key] " (" (features[key]/total)*100 "%)" }; print "Total: " total}' gencode.v31.annotation.gtf 
    # Feature exon: 1363843 (47.3247%)
    # Feature CDS: 755320 (26.2092%)
    # Feature UTR: 308315 (10.6984%)
    # Feature gene: 60603 (2.10289%)
    # Feature start_codon: 87299 (3.02923%)
    # Feature Selenocysteine: 119 (0.00412924%)
    # Feature stop_codon: 79505 (2.75878%)
    # Feature transcript: 226882 (7.87269%)
    # Total: 2881886
    ```

* Commands that were used in the comparison of `gene_types` for `gene` features between `GENCODE v21` and `GENCODE v31`:

    ```bash
    gawk '$3 ~ /gene/' gencode.v21.annotation.gtf | gawk 'match($0, /gene_type "([A-Za-z0-9]+)"/, a) { types[a[1]] += 1; total += 1 } END {for (key in types) { print "Gene type " key ": " types[key] " (" (types[key]/total)*100 "%)" }; print "Total: " total}'
    # Gene type rRNA: 549 (2.54508%)
    # Gene type antisense: 5542 (25.6919%)
    # Gene type pseudogene: 29 (0.13444%)
    # Gene type lincRNA: 7666 (35.5385%)
    # Gene type TEC: 1058 (4.90473%)
    # Gene type miRNA: 3837 (17.7878%)
    # Gene type snoRNA: 978 (4.53386%)
    # Gene type snRNA: 1912 (8.86375%)
    # Total: 21571
    
    gawk '$3 ~ /gene/' gencode.v31.annotation.gtf | gawk 'match($0, /gene_type "([A-Za-z0-9]+)"/, a) { types[a[1]] += 1; total += 1 } END {for (key in types) { print "Gene type " key ": " types[key] " (" (types[key]/total)*100 "%)" }; print "Total: " total}'
    # Gene type ribozyme: 8 (0.0351463%)
    # Gene type scaRNA: 49 (0.215271%)
    # Gene type lncRNA: 16840 (73.983%)
    # Gene type rRNA: 52 (0.228451%)
    # Gene type vaultRNA: 1 (0.00439329%)
    # Gene type scRNA: 1 (0.00439329%)
    # Gene type sRNA: 5 (0.0219664%)
    # Gene type pseudogene: 18 (0.0790792%)
    # Gene type TEC: 1064 (4.67446%)
    # Gene type miRNA: 1881 (8.26377%)
    # Gene type snoRNA: 942 (4.13848%)
    # Gene type snRNA: 1901 (8.35164%)
    # Total: 22762
    ```

* The following is a quick command that can be used on a GENCODE GTF to summarize the level counts/percentages in the GTF file. This was used to quantify how much of the GTF would be eliminated when all level 3 features were removed.

    ```bash
    gawk 'match($0, /level ([0-9]+)/, a) { levels[a[1]] += 1; total += 1 } END {for (key in levels) { print "Level " key ": " levels[key] " (" (levels[key]/total)*100 "%)" }; print "Total  : " total}' gencode.v31.annotation.gtf
    # Level 1: 186461 (6.4701%)
    # Level 2: 2450972 (85.0475%)
    # Level 3: 244453 (8.4824%)
    # Total  : 2881886
    ```
