<p align="center">
  <a href="https://github.com/stjudecloud/rfcs"><img src="https://github.com/stjudecloud/rfcs/raw/master/docs/rfcs-banner-blueprint.jpg" width="800" title="St. Jude Cloud RFCs"></a>
  <a href="https://travis-ci.org/stjudecloud/rfcs" target="_blank">
    <img alt="Build Status: Master" src="https://travis-ci.org/stjudecloud/rfcs.svg?branch=master" />
  </a>
  <a href="https://gitter.im/stjudecloud/rfcs?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge" target="_blank">
    <img alt="Gitter: RFCs" src="https://badges.gitter.im/stjudecloud/rfcs.svg" />
  </a>
</p>

> This repository contains all Request for Comments (or [rfcs][rfcs]) for the St. Jude Cloud project. Currently, RFCs on the St. Jude Cloud project are focused on outlining changes to the genomic analysis pipelines. However, we may expand to other areas over time.

### 🏠 [Homepage](https://stjudecloud.github.io/rfcs/)

## Process and Evaluation

Initiation of an RFC for the St. Jude Cloud project is currently limited to members of the St. Jude Cloud team. If you are a member of our community and have an idea you would like us to consider (or if you would like to propose your own RFC), please [contact us][contact] before investing a large amount of time in your proposal.

### Overview

All RFCs go through a three-stage process from draft to adoption:

1. An initial draft is constructed by the author and discussed internally with members of the core St. Jude Cloud team. 
   * Titles of pull requests in this phase are typical prefixed with `[WIP]` and should not be considered mature enough to accept comments from the community. 
   * The amount of time the RFC spends in this state is variable depending on the scope of the proposal.
2. The RFC opens for discussion from users and stakeholders within the St. Jude community. 
   * The beginning of this phase will be indicated by a comment on the pull request from the author.
   * The RFC will remain in this phase for **1 week**.
3. The RFC opens for discussion from the broader community.
   * The beginning of this phase will be indicated by a comment on the pull request from the author.
   * The RFC will remain in this phase for **2 weeks**.

Each phase of the process brings about further refinement of details and hardening of the proposal.

### Evaluation and acceptance

RFCs must have majority support within the St. Jude Cloud team before being adopted. In practice, most RFCs already have a significant amount of internal support before being drafted, so the rejection of an RFC will be exceedingly rare. The intent of these discussions is mostly geared towards gaining feedback from the community. 

We welcome all ideas and concerns as we seek to understand how the needs of the genomics community evolve over time. Our ultimate goal is to make informed decisions that will benefit the majority of our users. That said, the project's direction is not driven strictly by community consensus and ultimately will be decided by the St. Jude Cloud team based on our best judgement. 


### Creating a proposal

1. Create an initial draft of the RFC by creating a pull request on this repo.
   * Create a branch on this repo with the form `<username>/<featurename>`.
   * On that branch, copy the `0000-template.md` file to the `text/` folder and replace "template" with an appropriate feature name (e.g. `text/rnaseq-workflow-v2`).
   * If necessary, you can create a folder at `resources/<filename.md>/` to hold any images or supporting materials (e.g. `resources/0000-rnaseq-workflow-v2.0/`) 
   * Ensure that your pull request has `[WIP]` prefixing the title as long as it is not ready for review!
2. Open up the discussion for the community.
   * Ensure your pull request's main comment has a "rendered" tag for easy review (see [this example](https://github.com/stjudecloud/rfcs/pull/1)).
   * Remove the `[WIP]` if necessary.
   * Assign yourself as the assignee and add any relevant community members you'd like to review.
3. Once the RFC has reached a consensus, it is ready for inclusion in the main repository.
   * An owner of this repository will comment on the pull request with an assigned RFC number.
   * Change your filename and resources folder to reflect this change (e.g. rename `text/0000-rnaseq-workflow-v2.0.md` to `text/1234-rnaseq-workflow-v2.0.md`).
   * Ensure all of your links still work in your rendered copy.
   * Merge in the PR and delete the branch.

## Install

```sh
cargo install mdbook
```

## Usage

```sh
mdbook build
python3 -m http.server -d book
# visit the rendered version in your browser at http://localhost:8000.
```

## Author

👤 **St. Jude Cloud Team**

* Website: https://stjude.cloud
* Github: [@stjudecloud](https://github.com/stjudecloud)
* Twitter: [@StJudeResearch](https://twitter.com/StJudeResearch)

## Contributing

Contributions, issues and feature requests are welcome!<br />Feel free to check [issues page](https://github.com/stjudecloud/rfcs/issues). You can also take a look at the [contributing guide](https://github.com/stjudecloud/rfcs/blob/master/CONTRIBUTING.md).


## 📝 License

Copyright © 2020 [St. Jude Cloud Team](https://github.com/stjudecloud).<br />
This project is [MIT](https://github.com/stjudecloud/rfcs/blob/master/LICENSE.md) licensed.

## Questions

With any quetions, please [contact us][contact].

[rfcs]: https://en.wikipedia.org/wiki/Request_for_Comments
[contact]: mailto:support@stjude.cloud

