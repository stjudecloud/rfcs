# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
- [Specification](#specification)
- [Appendix](#appendix)

# Introduction

This RFC lays out the specification for the RNA-Seq mapping pipeline v2.0.0. The improvements contained within are largely based on (a) new version of tools/reference files and (b) feedback from the community. You can find the relevant discussion on the [associated pull request](https://github.com/stjudecloud/rfcs/pull/1).

# Motivation

- **Tool additions and updates.** The tools used in the first version of the RNA-seq pipeline are now woefully out of date (~2 years old). We desire to upgrade the tools to their latest published versions available so we can enjoy the benefits they provide. Furthermore, we'd like to explore adding many more tools to begin publishing QC results (probably a different QC pipeline in the future). See the [section below on tool additions and updates](#tool-additions-updates) for more details.
- **Updated reference files.** No changes have really been made to the `GRCh38_no_alt` analysis set reference genome. However, three major releases of the GENCODE gene model have transpired since we released the first revision of the RNA-Seq workflow. [GENCODE v31](https://www.gencodegenes.org/human/release_31.html) is the latest at the time of writing, so we'd like to utilize the new information contained there as well.
- **Quality of life improvements based on feedback from the community.** Input from the community has greatly informed our thoughts on what can be improved in this release. Some ideas that come to mind are:
  - Removal of the `ERCC SpikeIn` sequences to ensure consistent sequence dictionaries between DNA and RNA data.
    - Popular tools such as `GATK` and `picard` are generally unhappy if the sequence dictionaries don't match perfectly.
    - The inclusion of the External RNA Controls Consortium (ERCC) Spike-in Control RNA sequences in the alignment reference file we used for RNA-Seq mapping was hence causing issues when using mapped RNA-Seq BAM files in conjunction with other non-RNA-Seq BAM files in downstream analysis using these tools.
    - Last, many of our RNA-Seq samples were not generated using 'ERCC' spike-in control sequences, so the benefit its inclusion provides is minimal.

# Discussion

## Tool additions and upgrades

As part of the RNA-Seq workflow v2.0.0, multiple tools will be added and upgraded:

- `fq v0.2.0` ([Released](https://github.com/stjude/fqlib/releases/tag/v0.2.0) November 28, 2018) will be added. This tool will be used to validate the output of `picard SamToFastq`. `picard SamToFastq` does not currently catch all of the errors we wish to catch at this stage (such as duplicate read names in the FastQ file). Thus, we will leverage this tool to independently validate that the data is well-formed by our definition of that phrase.
- `ngsderive v1.0.2` ([Released](https://github.com/claymcleod/ngsderive/releases/tag/v1.0.2) March 3, 2020) will be added. This tool will be used to empirically infer strandedness of RNA-seq experiments.
- Update `STAR 2.5.3a` ([Released](https://github.com/alexdobin/STAR/releases/tag/2.5.3a) March 17, 2017) to `STAR 2.7.1a` ([Released](https://github.com/alexdobin/STAR/releases/tag/2.7.1a) May 15, 2019). Upgraded to receive the benefits of bug fixes and software optimizations.
- Update `samtools 1.4.0` ([Released](https://github.com/samtools/samtools/releases/tag/1.4) March 13, 2017) to `samtools 1.9` ([Released](https://github.com/samtools/samtools/releases/tag/1.9) July 18, 2018). Updating the samtools version whenever possible is of particular interest to me due to the historical fragility of the samtools code (although it has seemed to get better over the last year or so).
- Update `picard 2.9.4` ([Released](https://github.com/broadinstitute/picard/releases/tag/2.9.4) June 15, 2017) to `picard 2.20.2` ([Released](https://github.com/broadinstitute/picard/releases/tag/2.20.2) May 28, 2019). Upgraded to receive the benefits of bug fixes and software optimizations.

## GENCODE compatability

One of the major discussions during this round of revisions was the compatability of the `GRCh38_no_alt` reference genome with the latest `GENCODE` gene set. It was posed as the following question:

> This has just been an outstanding question of mine for a while — how big of an impact (if any at all) does the mismatch between the patch builds have? GENCODE v31 is built against `GRCh38.p12` (see the in the title on the [webpage](https://www.gencodegenes.org/human/release_31.html)), but obviously the no alt analysis set is derived from `GRCh38.p0` (with the pseudoautosomal regions hard masked, the EBV chromosome added, etc.)

As is apparent in the question, the `GRCh38` reference genome is regularly patched with non-coordinate altering, backwards compatable changes (see more information [here](https://www.ncbi.nlm.nih.gov/grc/help/patches/)). On the other hand, each `GENCODE` gene model release is based on a particular patch of the reference genome. This may be problematic because most bioinformatics analyses use a reference genome based on the non-patched version of `GRCh38`. What follows is our discussion and investigation into the effect between mismatching nucleotide sequences in the reference genome and the gene model.

First, we researched what some of the projects we respect in the community are doing:

| Pipeline                                                                 | Reference Genome                                                     | Reference Genome Patch | Gene Model                 | Gene Model Patch |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------- | ---------------------- | -------------------------- | ---------------- |
| GDC's [mRNA-Seq pipeline][gdc-mrnaseq-pipeline]                          | [`GRCh38_no_alt`-based w/ decoys + viral][gdc-reference-genome]      | `GRCh38.p0`            | [GENCODE v22][gencode-v22] | `GRCh38.p2`      |
| ENCODE's [RNA-Seq pipeline][encode-rnaseq-pipeline]                      | [`GRCh38_no_alt`-based w/ SpikeIns][encode-reference-genome]         | `GRCh38.p0`            | [GENCODE v24][gencode-v24] | `GRCh38.p5`      |
| Broad Institute's [GTEx + TOPMed RNA-Seq pipeline][gtex-rnaseq-pipeline] | [Broad's `GRCh38` w/ ERCC SpikeIn][broad-institute-reference-genome] | `GRCh38.p0`            | [GENCODE v26][gencode-v26] | `GRCh38.p10`     |

[gdc-mrnaseq-pipeline]: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
[gdc-reference-genome]: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
[encode-rnaseq-pipeline]: https://www.encodeproject.org/pipelines/ENCPL002LPE/https://www.encodeproject.org/pages/pipelines/#RNA-seq
[encode-reference-genome]: https://www.encodeproject.org/files/ENCFF742NER/
[gtex-rnaseq-pipeline]: https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq#reference-genome-and-annotation
[broad-institute-reference-genome]: https://software.broadinstitute.org/gatk/download/bundle
[gencode-v22]: https://www.gencodegenes.org/human/release_22.html
[gencode-v24]: https://www.gencodegenes.org/human/release_24.html
[gencode-v26]: https://www.gencodegenes.org/human/release_26.html

**Note:** You can confirm which patch the GENCODE genesets is based on just by clicking on the hyperlink. Verifying that each of these reference genomes is really based on `GRCh38_no_alt` takes a little bit more elbow grease: if you're interested, you can check out the comparison table [in the appendix](#reference-genome-comparison). If you are _really_ interested, you can recapitulate those results by running [the associated Jupyter notebook](https://github.com/stjudecloud/rfcs/tree/master/resources/0001-rnaseq-workflow-v2.0.0/GenomeComparison.ipynb).

Based on the results of the above investigation, I reached out to the author of STAR, Alex Dobin, to get his opinion on whether the differences might affect some results. You can read my question and his reply [here](https://github.com/alexdobin/STAR/issues/673). In short, he confirms that, yes, this may alter results for the various analysis types we were interested in.

Given the landscape of the community and the author's response, we considered three possible options for moving ahead:

1. We could try using the reference FASTA supplied with the respective GENCODE release as suggested by Dr. Dobin. This was the most undesirable approach from our groups perspective:
   - The main reason was that this reference genome did not mirror the sequence dictionary or characteristics of the reference genome we use in our DNA-Seq pipeline out of the box. As outlined in [the motivation section](#Motivation) of this document, this incongruency was a major driving factor for the creation of this RFC.
   - We agreed it would require too large an amount of postprocessing of the GENCODE reference genome to convert to apply all of the no alt analysis set changes (e.g. converting to UCSC names, masking regions of the genome, inserting the EBV chromosome).
   - Additionally, this could leave room for the introduction of lots of strange errors, and there was no interest in getting into the business of genome generation (there is a reason no one applies the patches to their no alt analysis set).
2. We could downgrade the GENCODE gene model to the latest release that was still based on `GRCh38.p0` ([GENCODE v21](https://www.gencodegenes.org/human/release_21.html) would the correct version to use in this case).
   - The concordance of the reference sequences obviously made this choice an attractive option. However, `GENCODE v21` was released over 5 years ago (06.2014) and there were many valuable updates over that time. In particular, a quick look showed that there were many more transcripts and that the number of lncRNAs more than doubled (see appendix). We did not want to lose all of this forward progress if we could avoid it. To see what we looked at, you can see [the relevant section in the appendix](#GENCODE-feature-comparisons).
3. We could use the latest GENCODE release despite the mismatches outlined above as it appears most other projects have done.
   - The general consensus was that the effect of this discordance would be small enough to tolerate so that we could gain all of the knowledge accumulated since the older release.
   - To quantify this, we measured the differences in gene expression (as measured by the `R^2` value between `GENCODE v21` and `GENCODE v31`) and splice junction detection (as measured by the number of splice junctions detected and the relative proportions of novel/partially novel/known junctions).

After discussion with the group, we decided to stick with option #3.

## Gene model post-processing

Originally, I had posed this question to the group:

> - Previously, we were filtering out anything not matching "level 1" or "level 2" from the gene model. This was due to best practices outlined during our RNA-Seq Workflow v1.0 discussions. I propose we revert this for the following reasons:
>   - The first sentence in section 2.2.2 of the [STAR 2.7.1.a manual](https://github.com/alexdobin/STAR/blob/2.7.1a/doc/STARmanual.pdf): "The use of the most comprehensive annotations for a given species is strongly recommended". So it seems the author recommends you use the most comprehensive gene model.
>   - Here is what [the GENCODE FAQ](https://www.gencodegenes.org/pages/faq.html) has to say about the level 3 annotations: "Ensembl loci where they are different from the Havana annotation or where no Havana annotation can be found". Given that the GENCODE geneset is the union of automated annotations from the `Ensembl-genebuild` and manual curation of the `Ensembl-Havana` team, this level should be harmless in the event that levels 1 & 2 don't apply.
>   - Last, the various other pipelines in the community don't tend to remove these features:
>     - The [GDC documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#rna-seq-alignment-command-line-parameters) does not currently describe any filtering of their [GENCODE v22][gencode-v22] GTF (see the section labeled "Step 1: Building the STAR index").
>     - ENCODE is not filtering any features, they are just changing some of the contig names in the [GENCODE v24][gencode-v24] GTF (see the analysis [in the appendix](#ENCODE-GTF-generation)).
>     - The TOPMed pipeline [does outline some postprocessing](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md#reference-annotation) they are doing to the [GENCODE v26][gencode-v26] GTF, but it's mostly just collapsing exons. I inspected the script, which does have some functionality built in to blacklist a list of transcript IDs. However, the documentation does not specify that this is turned on by default.

After discussion internally, we decided to discontinue removing `level 3` annotations by default. This is more consistent with what is being done in the community and it was decided that this was the most straightforward method with little associated risk. Therefore we are no longer performing any post-processing of the gene model.

## Quality of life improvements

- For dependency management, we have moved to using `conda` within standard docker images. All packages should be available within the `defaults`, `conda-forge`, and `bioconda` repositories.
- Add a checksum algorithm and publish the results in the data browser. After consideration internally, we decided to use the `md5sum` checksum because of its ubiquity.

## Various other changes

- Removed a section of the pipeline that reformatted the header of the BAM to be cleaner. `STAR` outputs a header that is formatted fine already, and I found this code to just be an area where an error could be introduced for little benefit.
- Removed a section of custom code that checks for duplicate read groups. `picard ValidateSamFile` does this for you (see [the documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=7571) for this tool. Specifically, the `DUPLICATE_READ_GROUP_ID` error).

# Specification

## Dependencies

If you'd like the full `conda` environment, you can install it using the following command. Obviously, you'll need to install [anaconda](https://www.anaconda.com/) first.

```bash
conda create -n star-mapping \
    -c conda-forge \
    -c bioconda \
    picard==2.20.2 \
    samtools==1.9 \
    star==2.7.1a \
    htseq==0.11.2 \
    -y
```

Additionally, you will want to install our `fqlib` library to check that FastQ files are properly paired and have no duplicate read names. Installation of the [Rust](https://rustup.rs/) programming language is required.

```bash
cargo install --git https://github.com/stjude/fqlib.git --tag v0.3.1
```

## Reference files

The following reference files are used as the basis of the RNA-Seq Workflow v2.0.0:

- Similarly to all analysis pipelines in St. Jude Cloud, we use the `GRCh38_no_alt` analysis set for our reference genome. You can get a copy of the file [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). Additionally, you can get the file by running the following commands:

  ```bash
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -O GRCh38_no_alt.fa.gz
  gunzip GRCh38_no_alt.fa.gz

  echo "a6da8681616c05eb542f1d91606a7b2f  GRCh38_no_alt.fa" > GRCh38_no_alt.fa.md5
  md5sum -c GRCh38_no_alt.fa.md5
  # > GRCh38_no_alt.fa: OK
  ```

- For the gene model, we use the GENCODE v31 "comprehensive gene annotation" GTF for the "CHR" regions. You can get a copy of the gene annotation file [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz). For the exact steps to generate the gene model we use, you can run the following commands:

  ```bash
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
  echo "0aac86e8a6c9e5fa5afe2a286325da81  gencode.v31.annotation.gtf.gz" > gencode.v31.annotation.gtf.gz.md5
  md5sum -c gencode.v31.annotation.gtf.gz.md5
  # > gencode.v31.annotation.gtf.gz: OK

  gunzip gencode.v31.annotation.gtf.gz
  echo "4e22351ae216e72aa57cd6d6011960f8  gencode.v31.annotation.gtf" > gencode.v31.annotation.gtf.md5
  md5sum -c gencode.v31.annotation.gtf.md5
  # > gencode.v31.annotation.gtf: OK
  ```

- Last, the following command is used to prepare the STAR index file:

  ```bash
  STAR --runMode genomeGenerate \                    # Use genome generation runMode.
       --genomeDir $OUTPUT_DIR \                     # Specify an output directory.
       --runThreadN $NCPU \                          # Number of threads to use to build genome database.
       --genomeFastaFiles $FASTA \                   # A path to the GRCh38_no_alt.fa FASTA file.
       --sjdbGTFfile $GENCODE_GTF_V31 \              # GENCODE v31 gene model file.
       --sjdbOverhang 125                            # Splice junction database overhang parameter, the optimal value is (Max length of RNA-Seq read-1).
  ```

## Workflow

Here are the resulting steps in the RNA-Seq Workflow v2.0.0 pipeline. There might be slight alterations in the actual implementation, which can be found in [the St. Jude Cloud workflows repository](https://github.com/stjudecloud/workflows/blob/master/workflows/rnaseq/rnaseq-standard.wdl).

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

5. Run the `STAR` alignment algorithm.

   ```bash
   STAR --readFilesIn $ALL_FASTQ_R1 $ALL_FASTQ_READ2 \ # FastQ files, separated by comma if there are multiple. The order of your R1 and R2 files has to match!
        --outSAMattrRGline $ALL_RG \                   # Read group lines in the same order as `readFilesIn` (derived from earlier `samtools split` step).
        --genomeDir $STARDB \                          # Directory containing the STAR genome
        --runThreadN $NCPU \                           # Number of threads to use. You must request the correct amount from HPCF first!
        --outSAMunmapped Within \                      # Keep unmapped reads in the final BAM file.
        --outSAMstrandField intronMotif \              # Preserve compatibility with Cufflinks by including the XS attribute (strand derived from intron motif).
        --outSAMattributes NH HI AS nM NM MD XS \      # Recommended SAM attributes to include for compatibility. Refer to manual for specifics.
        --outFilterMultimapScoreRange 1 \              # Ensures that all multi-mapped reads will need to share the mapping score.
        --outFilterMultimapNmax 20 \                   # Max number of multi-mapped SAM entries for each read.
        --outFilterMismatchNmax 10 \                   # Max number of mismatches allowed from alignment.
        --alignIntronMax 500000 \                      # Max intron size considered.
        --alignMatesGapMax 1000000 \                   # Max gap allowed between two mates.
        --sjdbScore 2 \                                # Additional weight given to alignments that cross database junctions when scoring alignments.
        --alignSJDBoverhangMin 1 \                     # Minimum overhang required on each side of a splice junction. Here, we require only one base on each side.
        --outFilterMatchNminOverLread 0.66 \           # 66% of the read must be perfectly matched to the reference sequence.
        --outFilterScoreMinOverLread 0.66 \            # Score must be greater than 66% of the read length. So for RL=100, the alignment must have a score > 66.
        --limitBAMsortRAM $RAM_LIMIT \                 # Amount of RAM to use for sorting. Recommended value is [Max amount of RAM] - 5GB.
        --outFileNamePrefix $OUT_FILE_PREFIX \         # All output files will have this path prepended.
        --twopassMode Basic                            # Use STAR two-pass mapping technique (refer to manual).
   ```

6. Run `picard SortSam` on the `STAR`-aligned BAM file. Note that this is much more memory efficient than using `STAR`'s built-in sorting (which often takes 100GB+ of RAM).

   ```bash
   picard SortSam I=$STAR_BAM \                  # Input BAM.
                  O=$MARKED_BAM \                # Duplicate-marked output BAM.
                  SO="coordinate" \              # 
                  CREATE_INDEX=false \           # Explicitly do not create an index at this step, in case the default changes.
                  CREATE_MD5_FILE=false \        # Explicity do not create an md5 checksum at this step, in case the default changes.
                  COMPRESSION_LEVEL=5 \          # Explicitly set the compression level to 5, although, at the time of writing, this is the default.
                  VALIDATION_STRINGENCY=SILENT   # Turn of validation stringency for this step.
   ```

7. Index the coordinate-sorted BAM file.

   ```bash
   samtools index $STAR_SORTED_BAM # STAR-aligned, coordinate-sorted BAM.
   ```

8. Run `picard ValidateSamFile` on the aligned and marked BAM file.

   ```bash
   picard ValidateSamFile I=$STAR_SORTED_BAM \    # STAR-aligned, coordinate-sorted BAM.
                  IGNORE=INVALID_PLATFORM_VALUE \ # Validations to ignore.
                  IGNORE=MISSING_PLATFORM_VALUE
   ```

9. Run `ngsderive`'s `strandedness` subcommand to confirm that the lab's information on strandedness reflects what was is computed. Manually triage any discrepancies. This is particularly useful for historical samples. Additionally,
if the value for strandedness isn't known at run time, we can use the inferred value (if reported).

   ```bash
   ngsderive strandedness $STAR_SORTED_BAM \     # STAR-aligned, coordinate-sorted BAM.
                          -g $GENCODE_GTF_V31_GZ # GENCODE v31 GTF (gzipped)
   ```

10. Next, `htseq-count` is run for the final counts file to be delivered.

    ```bash
    htseq-count -f bam \                            # Specify input file as BAM.
               -r pos \                             # Specify the BAM is position sorted.
               -s $PROVIDED_OR_INFERRED \           # Strandedness as specified by the lab and confirmed by `ngsderive strandedness` above. Typically `reverse` for St. Jude Cloud data.
               -m union \                           # For reference, GDC uses "intersection-nonempty".
               -i gene_name \                       # I'd like to use the colloquial gene name here. For reference, GDC uses "gene_id" here. Needs input from reviewers.
               --secondary-alignments ignore \      # Elect to ignore secondary alignments. Needs input from reviewers.
               --supplementary-alignments ignore \  # Elect to ignore supplementary alignments. Needs input from reviewers.
               $STAR_SORTED_BAM \                   # STAR-aligned, coordinate-sorted BAM.
               $GENCODE_GTF_V31                     # GENCODE v31 GTF
    ```

11. Generate the `md5sum` for the BAM

   ```bash
   md5sum $STAR_SORTED_BAM # STAR-aligned, coordinate-sorted BAM.
   ```

# Appendix

### Reference genome comparison

Below are the results of an analysis of each pipeline's `GRCh38`-based reference genome. The steps are essentially:

1. Research what reference genome each project is using and where to find it.
2. Download it.
3. Use `picard CreateSequenceDictionary` to get the md5sums for each sequence in the dictionary.
4. Compare for the common chromosomes in each reference (the autosomes, the sex chromosomes, and the EBV decoy).

If you are interested in seeing the _full_ comparison table or in regenerating these results, you can see [the associated Jupyter notebook](https://github.com/stjudecloud/rfcs/tree/master/resources/0001-rnaseq-workflow-v2.0.0/GenomeComparison.ipynb).

| Sequence Name | NCBI (baseline)                    | ENCODE                             | GDC                                | TOPMed                             | Concordant |
| ------------- | ---------------------------------- | ---------------------------------- | ---------------------------------- | ---------------------------------- | ---------- |
| chr1          | `6aef897c3d6ff0c78aff06ac189178dd` | `6aef897c3d6ff0c78aff06ac189178dd` | `6aef897c3d6ff0c78aff06ac189178dd` | `6aef897c3d6ff0c78aff06ac189178dd` | True       |
| chr2          | `f98db672eb0993dcfdabafe2a882905c` | `f98db672eb0993dcfdabafe2a882905c` | `f98db672eb0993dcfdabafe2a882905c` | `f98db672eb0993dcfdabafe2a882905c` | True       |
| chr3          | `76635a41ea913a405ded820447d067b0` | `76635a41ea913a405ded820447d067b0` | `76635a41ea913a405ded820447d067b0` | `76635a41ea913a405ded820447d067b0` | True       |
| chr4          | `3210fecf1eb92d5489da4346b3fddc6e` | `3210fecf1eb92d5489da4346b3fddc6e` | `3210fecf1eb92d5489da4346b3fddc6e` | `3210fecf1eb92d5489da4346b3fddc6e` | True       |
| chr5          | `a811b3dc9fe66af729dc0dddf7fa4f13` | `a811b3dc9fe66af729dc0dddf7fa4f13` | `a811b3dc9fe66af729dc0dddf7fa4f13` | `a811b3dc9fe66af729dc0dddf7fa4f13` | True       |
| chr6          | `5691468a67c7e7a7b5f2a3a683792c29` | `5691468a67c7e7a7b5f2a3a683792c29` | `5691468a67c7e7a7b5f2a3a683792c29` | `5691468a67c7e7a7b5f2a3a683792c29` | True       |
| chr7          | `cc044cc2256a1141212660fb07b6171e` | `cc044cc2256a1141212660fb07b6171e` | `cc044cc2256a1141212660fb07b6171e` | `cc044cc2256a1141212660fb07b6171e` | True       |
| chr8          | `c67955b5f7815a9a1edfaa15893d3616` | `c67955b5f7815a9a1edfaa15893d3616` | `c67955b5f7815a9a1edfaa15893d3616` | `c67955b5f7815a9a1edfaa15893d3616` | True       |
| chr9          | `6c198acf68b5af7b9d676dfdd531b5de` | `6c198acf68b5af7b9d676dfdd531b5de` | `6c198acf68b5af7b9d676dfdd531b5de` | `6c198acf68b5af7b9d676dfdd531b5de` | True       |
| chr10         | `c0eeee7acfdaf31b770a509bdaa6e51a` | `c0eeee7acfdaf31b770a509bdaa6e51a` | `c0eeee7acfdaf31b770a509bdaa6e51a` | `c0eeee7acfdaf31b770a509bdaa6e51a` | True       |
| chr11         | `1511375dc2dd1b633af8cf439ae90cec` | `1511375dc2dd1b633af8cf439ae90cec` | `1511375dc2dd1b633af8cf439ae90cec` | `1511375dc2dd1b633af8cf439ae90cec` | True       |
| chr12         | `96e414eace405d8c27a6d35ba19df56f` | `96e414eace405d8c27a6d35ba19df56f` | `96e414eace405d8c27a6d35ba19df56f` | `96e414eace405d8c27a6d35ba19df56f` | True       |
| chr13         | `a5437debe2ef9c9ef8f3ea2874ae1d82` | `a5437debe2ef9c9ef8f3ea2874ae1d82` | `a5437debe2ef9c9ef8f3ea2874ae1d82` | `a5437debe2ef9c9ef8f3ea2874ae1d82` | True       |
| chr14         | `e0f0eecc3bcab6178c62b6211565c807` | `e0f0eecc3bcab6178c62b6211565c807` | `e0f0eecc3bcab6178c62b6211565c807` | `e0f0eecc3bcab6178c62b6211565c807` | True       |
| chr15         | `f036bd11158407596ca6bf3581454706` | `f036bd11158407596ca6bf3581454706` | `f036bd11158407596ca6bf3581454706` | `f036bd11158407596ca6bf3581454706` | True       |
| chr16         | `db2d37c8b7d019caaf2dd64ba3a6f33a` | `db2d37c8b7d019caaf2dd64ba3a6f33a` | `db2d37c8b7d019caaf2dd64ba3a6f33a` | `db2d37c8b7d019caaf2dd64ba3a6f33a` | True       |
| chr17         | `f9a0fb01553adb183568e3eb9d8626db` | `f9a0fb01553adb183568e3eb9d8626db` | `f9a0fb01553adb183568e3eb9d8626db` | `f9a0fb01553adb183568e3eb9d8626db` | True       |
| chr18         | `11eeaa801f6b0e2e36a1138616b8ee9a` | `11eeaa801f6b0e2e36a1138616b8ee9a` | `11eeaa801f6b0e2e36a1138616b8ee9a` | `11eeaa801f6b0e2e36a1138616b8ee9a` | True       |
| chr19         | `85f9f4fc152c58cb7913c06d6b98573a` | `85f9f4fc152c58cb7913c06d6b98573a` | `85f9f4fc152c58cb7913c06d6b98573a` | `85f9f4fc152c58cb7913c06d6b98573a` | True       |
| chr20         | `b18e6c531b0bd70e949a7fc20859cb01` | `b18e6c531b0bd70e949a7fc20859cb01` | `b18e6c531b0bd70e949a7fc20859cb01` | `b18e6c531b0bd70e949a7fc20859cb01` | True       |
| chr21         | `974dc7aec0b755b19f031418fdedf293` | `974dc7aec0b755b19f031418fdedf293` | `974dc7aec0b755b19f031418fdedf293` | `974dc7aec0b755b19f031418fdedf293` | True       |
| chr22         | `ac37ec46683600f808cdd41eac1d55cd` | `ac37ec46683600f808cdd41eac1d55cd` | `ac37ec46683600f808cdd41eac1d55cd` | `ac37ec46683600f808cdd41eac1d55cd` | True       |
| chrX          | `2b3a55ff7f58eb308420c8a9b11cac50` | `2b3a55ff7f58eb308420c8a9b11cac50` | `2b3a55ff7f58eb308420c8a9b11cac50` | `2b3a55ff7f58eb308420c8a9b11cac50` | True       |
| chrY          | `ce3e31103314a704255f3cd90369ecce` | `ce3e31103314a704255f3cd90369ecce` | `ce3e31103314a704255f3cd90369ecce` | `ce3e31103314a704255f3cd90369ecce` | True       |
| chrM          | `c68f52674c9fb33aef52dcf399755519` | `c68f52674c9fb33aef52dcf399755519` | `c68f52674c9fb33aef52dcf399755519` | `c68f52674c9fb33aef52dcf399755519` | True       |
| chrEBV        | `6743bd63b3ff2b5b8985d8933c53290a` | `6743bd63b3ff2b5b8985d8933c53290a` | `6743bd63b3ff2b5b8985d8933c53290a` | `6743bd63b3ff2b5b8985d8933c53290a` | True       |

### ENCODE GTF generation

A small investigation was done to see how ENCODE's GENCODE v24 GTF was being produced (the md5sum between the original GENCODE v24 primary assembly file and ENCODE's version of the file weren't matching up). In short, they just changed the names of some of the contigs. You can verify this yourself by doing an `md5sum` without the sequence name column of the GTF (or manually download and inspect the differences like I did :)).

```bash
curl -sL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.primary_assembly.annotation.gtf.gz | \
         gunzip -c | \
         cut -f2- -d$'\t' | \
         md5sum
# 5227fac91c8d32b3fa3a8f78c4bf0e5c  -

curl -sL https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz | \
         gunzip -c | \
         cut -f2- -d$'\t' | \
         md5sum
# 5227fac91c8d32b3fa3a8f78c4bf0e5c  -
```

### GENCODE feature comparisons

Below are a few commands used to quickly evaluate how much the GENCODE geneset has changed over time. This was useful in our discussion about how much value would be lost if we just used `GENCODE v21` (which was based on `GRCh38.p0`). See [the discussion on the GENCODE compatibility](#GENCODE-compatability) for more information.

- Commands that were used in the comparison of feature types between `GENCODE v21` and `GENCODE v31`:

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
  # Total: 2546594

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

- Commands that were used in the comparison of `gene_types` for `gene` features between `GENCODE v21` and `GENCODE v31`:

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

- The following is a quick command that can be used on a GENCODE GTF to summarize the level counts/percentages in the GTF file. This was used to quantify how much of the GTF would be eliminated when all level 3 features were removed.

  ```bash
  gawk 'match($0, /level ([0-9]+)/, a) { levels[a[1]] += 1; total += 1 } END {for (key in levels) { print "Level " key ": " levels[key] " (" (levels[key]/total)*100 "%)" }; print "Total  : " total}' gencode.v31.annotation.gtf
  # Level 1: 186461 (6.4701%)
  # Level 2: 2450972 (85.0475%)
  # Level 3: 244453 (8.4824%)
  # Total  : 2881886
  ```
