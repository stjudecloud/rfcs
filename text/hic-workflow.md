# Hi-C Pipeline <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
- [Specification](#specification)

# Introduction

This RFC lays out the specification for the Hi-C harmonization pipeline. 

# Motivation

To provide the  community access to data from Hi-C experiments performed at St. Jude Children's Research Hospital, we propose the following data harmonization pipeline. The goal of this pipeline is to provide harmonized files for Hi-C data. For this pipeline, we will make no recommendations on downstream analysis, focusing instead on harmonizing the underlying sequencing data and leaving analysis decisions to the user.

# Discussion

Hi-C investigates chromatin conformation. 

# Specification

## Dependencies

- Picard Tools

## Reference Files

There are no reference files for the Hi-C harmonization workflow.


## Workflow

Here are the resulting steps in the Hi-C workflow. There might be slight alterations in the actual implementation, which can be found in [the St. Jude Cloud workflows repository](https://github.com/stjudecloud/workflows/blob/master/workflows/hic/hic-standard.wdl).

1. Run `picard ValidateSam` on the incoming BAM to ensure that it is well-formed enough to strip alignment information.

    ```bash
    picard ValidateSamFile \
                        I=$INPUT_BAM \                     # Input BAM.
                        IGNORE=INVALID_PLATFORM_VALUE \    # Validations to ignore.
                        IGNORE=MISSING_PLATFORM_VALUE
    ```

2. Run Picard `RevertSam` on the validated BAM file.

    ```bash
        picard RevertSam \
                INPUT=$INPUT_BAM \                        # Input BMA to revert
                OUTPUT=$REVERTED_BAM \                    # Output unaligned BAM name
                REMOVE_ALIGNMENT_INFORMATION=true \       # Remove alignments
                REMOVE_DUPLICATE_INFORMATION=true \       # Remove duplicate flags
                VALIDATION_STRINGENCY=SILENT \            # Ignore some validation warnings
                SORT_ORDER=queryname                      # Sort by queryname
    ```