# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Pre-Discussion](#pre-discussion)
- [Discussion](#discussion)

# Introduction

This RFC lays out the Data Versioning process for St. Jude Cloud.

# Motivation

To provide the community with harmonized data from experiments performed at St. Jude Children's Research Hospital, the St. Jude Cloud team has developed a datatype-specific harmonization workflows. Over time, it is necessary to incorporate bug fixes and other changes in those workflows and the underlying tools. In order to make these changes with minimal impact, we seek to define "data versions", within which, we have evidence that the workflow results are compatible and not "substantially" different. This RFC will lay out the methods by which we determine what constitutes a significant difference in the resulting data.

# Pre-Discussion

We have previously had short discussions regarding changes to our data versioning scheme. Currently we attach the pipeline version used to generate the data to the file's metadata. This has left us in a position where we are hesitant to increment pipelines to new versions or update to new versions of underlying tools. We generally agree that there are many pipeline-level changes that do not materially impact the output results. The goal of this RFC is to determine what changes we can be made without creating a new data version. To do so, we also need to generate metrics and data that we can use to justify to users that the results are not different.

There are a number of approaches that could be taken to do this. We could look at alignment information, in which case, we may want to use results from the QC pipeline to justify. Another option is to examine the results from an analysis perspective. If the variant calls are the same (or reasonably the same), that may be sufficient. These are just two limited examples of options we could use.

Currently we have the RNA-Seq workflow v2.0.0 in use for St. Jude Cloud. We have subsequently released improvements, including the integration of XenoCP. This is not in use for production work as we include the pipeline version with the data files in St. Jude Cloud. Creation of a data versioning standard will enable us to roll out bug fixes and new features on an appropriate timeline and provide confidence to ourselves and our users that data is consistent, hopefully, without requiring complete reruns of all data for any change.

We also have issues with the `msgen` pipeline that is used for WGS and WES. Microsoft has made changes to the pipeline and process throughout the time that we have been using it for production data releases. We have not incremented any versions in that time. So as it stands, we are trusting, without any verification, that any changes introduced are acceptable. 

A data versioning process should help us find a reasonable balance that doesn't stop us from making any changes, ever, but also we want to avoid incorporating any and all changes without any assurances, or insight in `msgen`'s case, that those changes are reasonable and not meaningfully impactful to the data results.

# Discussion

## Whole-Genome and Whole-Exome Sequencing

We propose to evaluate whole-genome (WGS) and whole-exome (WES) sequencing by running well characterized samples through each iteration of the analysis pipeline. The resulting variant calls (gVCFs) will be compared to existing high-quality variant calls. This comparison will be conducted using Illumina's `hap.py` [comparison tool](https://github.com/Illumina/hap.py) as [recommended](https://www.biorxiv.org/content/10.1101/270157v3) by the Global Alliance for Genomics and Health (GA4GH) Benchmarking Team. Specifically, we propose to run samples from the National Institute for Standards and Technology (NIST)'s Genome in a Bottle (GIAB) project. We will  perform analysis using samples HG002, HG003, HG004, HG005, HG006, and HG007 for WGS. For WES, we will use samples HG002, HG003, HG004, and HG005. The results from prior iterations of the pipeline will be supplied as the truth set. The confident call sets from GIAB will be provided as the gold standard dataset. The variant calls from the new workflow version will be treated as the query.

## RNA-Seq

We propose to evaluate RNA sequencing (RNA-Seq) by running well characterized samples through each iteration of the analysis pipeline. Specifically, we propose to run samples from the National Institute for Standards and Technology (NIST)'s Genome in a Bottle (GIAB) project. We will use RNA-Seq data from HG002, HG004, and HG005 samples. These three samples will be run through new iterations of the [St. Jude Cloud RNA-Seq harmonization workflow](https://stjudecloud.github.io/rfcs/0001-rnaseq-workflow-v2.0.0.html). The pipelines will output aligned BAMs and feature count files for each sample. We will then generate variant calls using the [GATK RNA-Seq short variant discovery best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-). We will use Illumina's hap.py comparison tool, along with the high confidence variant calls from GIAB to compare output from prior versions of the RNA-Seq pipeline to the proposed version.

### GATK RNA-Seq short variant discovery best practices workflow

The GATK best practices workflows are de facto standards for data processing in bioinformatics. Therefore we will reuse their existing workflow for generating variant calls from RNA-Seq data. The workflow accepts a STAR-aligned BAM file, such as those produced by our RNA-Seq alignment harmonization workflow, and produces a set of variant calls by running a number of GATK steps.

![GATK RNA-Seq short variant discovery workflow diagram](../resources/data_versioning/rnaseq-variant-workflow.png)

#### Workflow

1. Run Picard `MarkDuplicates`

    ```bash
        picard -Xmx~{java_heap_size}g MarkDuplicates \
            -I ~{bam} \
            -O ~{if create_bam then prefix + ".bam" else "/dev/null"} \
            --VALIDATION_STRINGENCY SILENT \
            --CREATE_INDEX ~{create_bam} \
            --CREATE_MD5_FILE ~{create_bam} \
            --COMPRESSION_LEVEL 5 \
            --METRICS_FILE ~{prefix}.metrics.txt
    ```

2. Run GATK `SplitNCigarReads`

    ```bash
        gatk \
            SplitNCigarReads \
            -R ~{fasta} \
            -I ~{bam} \
            -O ~{prefix}.bam \
            -OBM true
       # GATK is unreasonable and uses the plain ".bai" suffix.
       mv ~{prefix}.bai ~{prefix}.bam.bai
    ```

3. Run GATK Spark-enabled `BaseRecalibrator` (`BaseRecalibratorSpark`).

   ```bash
        gatk \
            --java-options "-Xms4000m" \
            BaseRecalibratorSpark \
            -R ~{fasta} \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.txt \
            -known-sites ~{dbSNP_vcf} \
            -known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
            --spark-master local[4]
   ```

4. Run GATK Spark-enabled `ApplyBQSR` (`ApplyBQSRSpark`).

    ```bash
            gatk \
            --java-options "-Xms3000m" \
            ApplyBQSRSpark \
            --spark-master local[4] \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.bam \
            --bqsr-recal-file ~{recalibration_report}
    ```

5. Run Picard `IntervalListTools` to generate list of regions to scatter over.

    ```bash
        picard -Xms1g \
            IntervalListTools \
            SCATTER_COUNT=~{scatter_count} \
            SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            UNIQUE=~{unique} \
            SORT=~{sort} \
            INPUT=~{interval_list} \
            OUTPUT=out
    ```

6. Run GATK Spark-enabled `HaplotypeCaller` (`HaplotypeCallerSpark`) on each region.

    ```bash
        gatk \
           --java-options "-Xms6000m" \
            HaplotypeCallerSpark \
            -R ~{fasta} \
            -I ~{bam} \
            -L ~{interval_list} \
            -O ~{prefix}.vcf.gz \
            ~{if use_soft_clipped_bases then "" else "--dont-use-soft-clipped-bases"} \
            --standard-min-confidence-threshold-for-calling ~{stand_call_conf} \
           --spark-master local[4]
    ```

7. Run Picard `MergeVcfs`.

    ```bash
        picard -Xms2000m \
            MergeVcfs \
            ~{sep(' ', prefix('--INPUT=', vcfs))} \
            --OUTPUT ~{output_vcf_name}
    ```

8. Run GATK `VariantFiltration`.

    ```bash
        gatk \
            VariantFiltration \
                --R ~{fasta} \
                --V ~{vcf} \
                --window ~{window} \
                --cluster ~{cluster} \
                 ~{sep(' ', prefix('--filter-name ', filter_name))} \
                 ~{sep(' ', prefix('--filter-expression ', squote(filter_expression)))} \
                -O ~{prefix}.vcf.gz
    ```

### Validation

After generation of RNA-Seq variants, we will run Illumina's `hap.py` comparison tool, as recommended by GA4GH, to evaluate the variant calls produced by the newest version of the pipeline compared with those produced by the prior version of the pipeline and compared against the high-quality variant call benchmark set published by GIAB.