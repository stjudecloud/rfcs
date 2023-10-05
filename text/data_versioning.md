# Table of Contents <!-- omit in toc -->

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Pre-Discussion](#pre-discussion)

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