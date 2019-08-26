# Request for Comments

[![Build Status](https://travis-ci.org/stjudecloud/rfcs.svg?branch=master)](https://travis-ci.org/stjudecloud/rfcs)

This repository contains all Request for Comment documents (or [rfcs][rfcs]) for the St. Jude Cloud project. Currently, the scope for what warrants an RFC on the St. Jude Cloud project is primarily focused on changes to the genomic analysis workflows. However, we may expand this scope to other areas over time as we work towards being as transparent as possible.

## Contributing

First, we should note that the initiator of an RFC for the St. Jude Cloud project is currently limited to members of the St. Jude Cloud team. If you have an idea you would like us to consider (or if you would like to propose your own RFC, please [contact us][contact] before investing a large amount of time in your proposal).

The process for RFCs is relatively immature at the current time and is subject to change. The current lifecycle of an RFC on this project is as follows:

1. Proposal of the RFC using a pull request on this repo.
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

## Build process

To build the repo and see the rendered copy of all RFCs, you can follow these steps:

```bash
bash bin/generate.sh
python3 -m http.server -d book
# visit the rendered version in your browser at http://localhost:8000.
```

## Questions

With any quetions, please [contact us][contact].

[rfcs]: https://en.wikipedia.org/wiki/Request_for_Comments
[contact]: mailto:support@stjude.cloud
