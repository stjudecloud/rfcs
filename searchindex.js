Object.assign(window.search, {"doc_urls":["introduction.html#request-for-comments","introduction.html#contributing","introduction.html#build-process","introduction.html#questions"],"index":{"documentStore":{"docInfo":{"0":{"body":39,"breadcrumbs":2,"title":2},"1":{"body":158,"breadcrumbs":1,"title":1},"2":{"body":20,"breadcrumbs":2,"title":2},"3":{"body":3,"breadcrumbs":1,"title":1}},"docs":{"0":{"body":"Build Status Join the chat at https://gitter.im/stjudecloud/rfcs This repository contains all Request for Comment documents (or rfcs ) for the St. Jude Cloud project. Currently, the scope for what warrants an RFC on the St. Jude Cloud project is primarily focused on changes to the genomic analysis workflows. However, we may expand this scope to other areas over time as we work towards being as transparent as possible.","breadcrumbs":"Request for Comments","id":"0","title":"Request for Comments"},"1":{"body":"First, we should note that the initiator of an RFC for the St. Jude Cloud project is currently limited to members of the St. Jude Cloud team. If you have an idea you would like us to consider (or if you would like to propose your own RFC, please contact us before investing a large amount of time in your proposal). The process for RFCs is relatively immature at the current time and is subject to change. The current lifecycle of an RFC on this project is as follows: Proposal of the RFC using a pull request on this repo. Create a branch on this repo with the form <username>/<featurename>. On that branch, copy the 0000-template.md file to the text/ folder and replace \"template\" with an appropriate feature name (e.g. text/rnaseq-workflow-v2). If necessary, you can create a folder at resources/<filename.md>/ to hold any images or supporting materials (e.g. resources/0000-rnaseq-workflow-v2.0/) Ensure that your pull request has [WIP] prefixing the title as long as it is not ready for review! Open up the discussion for the community. Ensure your pull request's main comment has a \"rendered\" tag for easy review (see this example ). Remove the [WIP] if necessary. Assign yourself as the assignee and add any relevant community members you'd like to review. Once the RFC has reached a consensus, it is ready for inclusion in the main repository. An owner of this repository will comment on the pull request with an assigned RFC number. Change your filename and resources folder to reflect this change (e.g. rename text/0000-rnaseq-workflow-v2.0.md to text/1234-rnaseq-workflow-v2.0.md). Ensure all of your links still work in your rendered copy. Merge in the PR and delete the branch.","breadcrumbs":"Contributing","id":"1","title":"Contributing"},"2":{"body":"To build the repo and see the rendered copy of all RFCs, you can follow these steps: bash bin/generate.sh\npython3 -m http.server -d book\n# visit the rendered version in your browser at http://localhost:8000.","breadcrumbs":"Build process","id":"2","title":"Build process"},"3":{"body":"With any quetions, please contact us .","breadcrumbs":"Questions","id":"3","title":"Questions"}},"length":4,"save":true},"fields":["title","body","breadcrumbs"],"index":{"body":{"root":{"0":{"0":{"0":{"0":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"a":{"d":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{},"m":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}}},"n":{"a":{"df":0,"docs":{},"l":{"df":0,"docs":{},"y":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"df":0,"docs":{}},"p":{"df":0,"docs":{},"p":{"df":0,"docs":{},"r":{"df":0,"docs":{},"o":{"df":0,"docs":{},"p":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}},"r":{"df":0,"docs":{},"e":{"a":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}},"s":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"g":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.4142135623730951}},"e":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}},"b":{"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"h":{"df":1,"docs":{"2":{"tf":1.0}}}}},"df":0,"docs":{},"e":{"df":1,"docs":{"0":{"tf":1.0}},"f":{"df":0,"docs":{},"o":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"df":0,"docs":{},"n":{"/":{"df":0,"docs":{},"g":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"e":{".":{"df":0,"docs":{},"s":{"df":0,"docs":{},"h":{"df":1,"docs":{"2":{"tf":1.0}}}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}}},"o":{"df":0,"docs":{},"o":{"df":0,"docs":{},"k":{"df":1,"docs":{"2":{"tf":1.0}}}}},"r":{"a":{"df":0,"docs":{},"n":{"c":{"df":0,"docs":{},"h":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}},"df":0,"docs":{}}},"df":0,"docs":{},"o":{"df":0,"docs":{},"w":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"2":{"tf":1.0}}}}}}}},"u":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"d":{"df":2,"docs":{"0":{"tf":1.0},"2":{"tf":1.4142135623730951}}},"df":0,"docs":{}}}}},"c":{"df":0,"docs":{},"h":{"a":{"df":0,"docs":{},"n":{"df":0,"docs":{},"g":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.7320508075688772}}}},"t":{"df":1,"docs":{"0":{"tf":1.0}}}},"df":0,"docs":{}},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"d":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}},"df":0,"docs":{}}}},"o":{"df":0,"docs":{},"m":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}}}},"u":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}}},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"u":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"t":{"a":{"c":{"df":0,"docs":{},"t":{"df":2,"docs":{"1":{"tf":1.0},"3":{"tf":1.0}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"n":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"b":{"df":0,"docs":{},"u":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{}}}}},"p":{"df":0,"docs":{},"i":{"df":2,"docs":{"1":{"tf":1.4142135623730951},"2":{"tf":1.0}}}}},"r":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}}},"u":{"df":0,"docs":{},"r":{"df":0,"docs":{},"r":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.7320508075688772}}}}}}}}},"d":{"df":1,"docs":{"2":{"tf":1.0}},"e":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"df":0,"docs":{},"s":{"c":{"df":0,"docs":{},"u":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}}},"o":{"c":{"df":0,"docs":{},"u":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"0":{"tf":1.0}}}}}}}},"df":0,"docs":{}}},"df":0,"docs":{},"e":{".":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}},"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}},"x":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"p":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{},"p":{"a":{"df":0,"docs":{},"n":{"d":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}}},"f":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"df":1,"docs":{"1":{"tf":1.0}},"n":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}},"r":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"o":{"c":{"df":0,"docs":{},"u":{"df":0,"docs":{},"s":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":0,"docs":{},"l":{"d":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}},"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"w":{"df":2,"docs":{"1":{"tf":1.0},"2":{"tf":1.0}}}}}},"r":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"g":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"o":{"df":0,"docs":{},"m":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"h":{"df":0,"docs":{},"o":{"df":0,"docs":{},"l":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"t":{"df":0,"docs":{},"t":{"df":0,"docs":{},"p":{".":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"v":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"2":{"tf":1.0}}}}}}}}},":":{"/":{"/":{"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"c":{"a":{"df":0,"docs":{},"l":{"df":0,"docs":{},"h":{"df":0,"docs":{},"o":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{":":{"8":{"0":{"0":{"0":{"df":1,"docs":{"2":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}},"df":0,"docs":{}}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{},"s":{":":{"/":{"/":{"df":0,"docs":{},"g":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":0,"docs":{},"t":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{".":{"df":0,"docs":{},"i":{"df":0,"docs":{},"m":{"/":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":0,"docs":{},"j":{"df":0,"docs":{},"u":{"d":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"d":{"/":{"df":0,"docs":{},"r":{"df":0,"docs":{},"f":{"c":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}}},"df":0,"docs":{}},"df":0,"docs":{}}}}},"df":0,"docs":{}}},"df":0,"docs":{}}}}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}}}}}},"i":{"d":{"df":0,"docs":{},"e":{"a":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{},"m":{"a":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{},"m":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}}},"n":{"c":{"df":0,"docs":{},"l":{"df":0,"docs":{},"u":{"df":0,"docs":{},"s":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}},"v":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}},"j":{"df":0,"docs":{},"o":{"df":0,"docs":{},"i":{"df":0,"docs":{},"n":{"df":1,"docs":{"0":{"tf":1.0}}}}},"u":{"d":{"df":0,"docs":{},"e":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}}},"l":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"f":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"y":{"c":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}},"df":0,"docs":{}}},"m":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}},"n":{"df":0,"docs":{},"k":{"df":1,"docs":{"1":{"tf":1.0}}}}},"o":{"df":0,"docs":{},"n":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"m":{"a":{"df":0,"docs":{},"i":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"t":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}}}},"df":1,"docs":{"2":{"tf":1.0}},"e":{"df":0,"docs":{},"m":{"b":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}},"r":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"n":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}}}}},"df":0,"docs":{}},"o":{"df":0,"docs":{},"t":{"df":0,"docs":{},"e":{"df":1,"docs":{"1":{"tf":1.0}}}}},"u":{"df":0,"docs":{},"m":{"b":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{}}}},"o":{"df":0,"docs":{},"n":{"c":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"p":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.0}}}}},"v":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"0":{"tf":1.0}}}}},"w":{"df":0,"docs":{},"n":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}}},"p":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"s":{"df":2,"docs":{"1":{"tf":1.0},"3":{"tf":1.0}}}},"df":0,"docs":{}}},"o":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"b":{"df":0,"docs":{},"l":{"df":1,"docs":{"0":{"tf":1.0}}}},"df":0,"docs":{}}}}},"r":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"f":{"df":0,"docs":{},"i":{"df":0,"docs":{},"x":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"df":0,"docs":{},"m":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"i":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"df":0,"docs":{}}},"o":{"c":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":2,"docs":{"1":{"tf":1.0},"2":{"tf":1.0}}}}}},"df":0,"docs":{},"j":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"t":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}}},"p":{"df":0,"docs":{},"o":{"df":0,"docs":{},"s":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}}},"u":{"df":0,"docs":{},"l":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":2.0}}}}},"y":{"df":0,"docs":{},"t":{"df":0,"docs":{},"h":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"3":{"df":1,"docs":{"2":{"tf":1.0}}},"df":0,"docs":{}}}}}}},"q":{"df":0,"docs":{},"u":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"3":{"tf":1.0}}}}}}},"t":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"3":{"tf":1.0}}}}}}}}},"r":{"df":0,"docs":{},"e":{"a":{"c":{"df":0,"docs":{},"h":{"df":1,"docs":{"1":{"tf":1.0}}}},"d":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}},"df":0,"docs":{},"f":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}},"l":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"v":{"df":1,"docs":{"1":{"tf":1.0}}}}},"m":{"df":0,"docs":{},"o":{"df":0,"docs":{},"v":{"df":1,"docs":{"1":{"tf":1.0}}}}},"n":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"d":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":2,"docs":{"1":{"tf":1.4142135623730951},"2":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}},"p":{"df":0,"docs":{},"l":{"a":{"c":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"o":{"df":2,"docs":{"1":{"tf":1.4142135623730951},"2":{"tf":1.0}},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":0,"docs":{},"o":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.4142135623730951}}}}}}}}}},"q":{"df":0,"docs":{},"u":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"'":{"df":1,"docs":{"1":{"tf":1.0}}},"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.7320508075688772}}}}}}},"s":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"c":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"s":{"/":{"0":{"0":{"0":{"0":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"<":{"df":0,"docs":{},"f":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{".":{"df":0,"docs":{},"m":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}},"v":{"df":0,"docs":{},"i":{"df":0,"docs":{},"e":{"df":0,"docs":{},"w":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}}},"f":{"c":{"df":3,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":2.6457513110645907},"2":{"tf":1.0}}},"df":0,"docs":{}},"n":{"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"q":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}},"df":0,"docs":{}}},"s":{"c":{"df":0,"docs":{},"o":{"df":0,"docs":{},"p":{"df":0,"docs":{},"e":{"df":1,"docs":{"0":{"tf":1.4142135623730951}}}}}},"df":0,"docs":{},"e":{"df":0,"docs":{},"e":{"df":2,"docs":{"1":{"tf":1.0},"2":{"tf":1.0}}}},"t":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}},"e":{"df":0,"docs":{},"p":{"df":1,"docs":{"2":{"tf":1.0}}}},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"u":{"b":{"df":0,"docs":{},"j":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}},"df":0,"docs":{},"p":{"df":0,"docs":{},"p":{"df":0,"docs":{},"o":{"df":0,"docs":{},"r":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}},"t":{"a":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{},"m":{"df":0,"docs":{},"p":{"df":0,"docs":{},"l":{"a":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}},"e":{".":{"df":0,"docs":{},"m":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}},"x":{"df":0,"docs":{},"t":{"/":{"0":{"0":{"0":{"0":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"1":{"2":{"3":{"4":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{},"r":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"q":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}}}},"df":1,"docs":{"1":{"tf":1.0}}}}},"i":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.4142135623730951}}}},"t":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}}},"o":{"df":0,"docs":{},"w":{"a":{"df":0,"docs":{},"r":{"d":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}},"r":{"a":{"df":0,"docs":{},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"p":{"a":{"df":0,"docs":{},"r":{"df":1,"docs":{"0":{"tf":1.0}}}},"df":0,"docs":{}}}}},"df":0,"docs":{}}},"u":{"df":0,"docs":{},"p":{"df":1,"docs":{"1":{"tf":1.0}}},"s":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{">":{"/":{"<":{"df":0,"docs":{},"f":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}},"v":{"2":{".":{"0":{".":{"df":0,"docs":{},"m":{"d":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}},"df":0,"docs":{}}},"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"2":{"tf":1.0}}}}}}}},"i":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":1,"docs":{"2":{"tf":1.0}}}}}}},"w":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"r":{"a":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":0,"docs":{}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"p":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"o":{"df":0,"docs":{},"r":{"df":0,"docs":{},"k":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.0}},"f":{"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"w":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":2.0}}}}}}}}}},"y":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"'":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{},"r":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"l":{"df":0,"docs":{},"f":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}}}}},"breadcrumbs":{"root":{"0":{"0":{"0":{"0":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"a":{"d":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{},"m":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}}},"n":{"a":{"df":0,"docs":{},"l":{"df":0,"docs":{},"y":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"df":0,"docs":{}},"p":{"df":0,"docs":{},"p":{"df":0,"docs":{},"r":{"df":0,"docs":{},"o":{"df":0,"docs":{},"p":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}},"r":{"df":0,"docs":{},"e":{"a":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}},"s":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"g":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.4142135623730951}},"e":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}},"b":{"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"h":{"df":1,"docs":{"2":{"tf":1.0}}}}},"df":0,"docs":{},"e":{"df":1,"docs":{"0":{"tf":1.0}},"f":{"df":0,"docs":{},"o":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"df":0,"docs":{},"n":{"/":{"df":0,"docs":{},"g":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"e":{".":{"df":0,"docs":{},"s":{"df":0,"docs":{},"h":{"df":1,"docs":{"2":{"tf":1.0}}}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}}},"o":{"df":0,"docs":{},"o":{"df":0,"docs":{},"k":{"df":1,"docs":{"2":{"tf":1.0}}}}},"r":{"a":{"df":0,"docs":{},"n":{"c":{"df":0,"docs":{},"h":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}},"df":0,"docs":{}}},"df":0,"docs":{},"o":{"df":0,"docs":{},"w":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"2":{"tf":1.0}}}}}}}},"u":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"d":{"df":2,"docs":{"0":{"tf":1.0},"2":{"tf":1.7320508075688772}}},"df":0,"docs":{}}}}},"c":{"df":0,"docs":{},"h":{"a":{"df":0,"docs":{},"n":{"df":0,"docs":{},"g":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.7320508075688772}}}},"t":{"df":1,"docs":{"0":{"tf":1.0}}}},"df":0,"docs":{}},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"d":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}},"df":0,"docs":{}}}},"o":{"df":0,"docs":{},"m":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":2,"docs":{"0":{"tf":1.7320508075688772},"1":{"tf":1.4142135623730951}}}}},"u":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}}},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"u":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"t":{"a":{"c":{"df":0,"docs":{},"t":{"df":2,"docs":{"1":{"tf":1.0},"3":{"tf":1.0}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"n":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"b":{"df":0,"docs":{},"u":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}}}}},"p":{"df":0,"docs":{},"i":{"df":2,"docs":{"1":{"tf":1.4142135623730951},"2":{"tf":1.0}}}}},"r":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}}},"u":{"df":0,"docs":{},"r":{"df":0,"docs":{},"r":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.7320508075688772}}}}}}}}},"d":{"df":1,"docs":{"2":{"tf":1.0}},"e":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"df":0,"docs":{},"s":{"c":{"df":0,"docs":{},"u":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}}},"o":{"c":{"df":0,"docs":{},"u":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"0":{"tf":1.0}}}}}}}},"df":0,"docs":{}}},"df":0,"docs":{},"e":{".":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}},"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}},"x":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"p":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{},"p":{"a":{"df":0,"docs":{},"n":{"d":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}}},"f":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"df":1,"docs":{"1":{"tf":1.0}},"n":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}},"r":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"o":{"c":{"df":0,"docs":{},"u":{"df":0,"docs":{},"s":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":0,"docs":{},"l":{"d":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}},"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"w":{"df":2,"docs":{"1":{"tf":1.0},"2":{"tf":1.0}}}}}},"r":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"g":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"o":{"df":0,"docs":{},"m":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"h":{"df":0,"docs":{},"o":{"df":0,"docs":{},"l":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"t":{"df":0,"docs":{},"t":{"df":0,"docs":{},"p":{".":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"v":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"2":{"tf":1.0}}}}}}}}},":":{"/":{"/":{"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"c":{"a":{"df":0,"docs":{},"l":{"df":0,"docs":{},"h":{"df":0,"docs":{},"o":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{":":{"8":{"0":{"0":{"0":{"df":1,"docs":{"2":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}},"df":0,"docs":{}}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{},"s":{":":{"/":{"/":{"df":0,"docs":{},"g":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":0,"docs":{},"t":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{".":{"df":0,"docs":{},"i":{"df":0,"docs":{},"m":{"/":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":0,"docs":{},"j":{"df":0,"docs":{},"u":{"d":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"d":{"/":{"df":0,"docs":{},"r":{"df":0,"docs":{},"f":{"c":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}}},"df":0,"docs":{}},"df":0,"docs":{}}}}},"df":0,"docs":{}}},"df":0,"docs":{}}}}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}}}}}},"i":{"d":{"df":0,"docs":{},"e":{"a":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{},"m":{"a":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{},"m":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}}},"n":{"c":{"df":0,"docs":{},"l":{"df":0,"docs":{},"u":{"df":0,"docs":{},"s":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}},"v":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}},"j":{"df":0,"docs":{},"o":{"df":0,"docs":{},"i":{"df":0,"docs":{},"n":{"df":1,"docs":{"0":{"tf":1.0}}}}},"u":{"d":{"df":0,"docs":{},"e":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}}},"l":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"f":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"y":{"c":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}},"df":0,"docs":{}}},"m":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}},"n":{"df":0,"docs":{},"k":{"df":1,"docs":{"1":{"tf":1.0}}}}},"o":{"df":0,"docs":{},"n":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"m":{"a":{"df":0,"docs":{},"i":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"t":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.0}}}}}}},"df":1,"docs":{"2":{"tf":1.0}},"e":{"df":0,"docs":{},"m":{"b":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}},"r":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"n":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}}}}},"df":0,"docs":{}},"o":{"df":0,"docs":{},"t":{"df":0,"docs":{},"e":{"df":1,"docs":{"1":{"tf":1.0}}}}},"u":{"df":0,"docs":{},"m":{"b":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{}}}},"o":{"df":0,"docs":{},"n":{"c":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"p":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":1,"docs":{"1":{"tf":1.0}}}}},"v":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"0":{"tf":1.0}}}}},"w":{"df":0,"docs":{},"n":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":1,"docs":{"1":{"tf":1.0}}}}}}},"p":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"s":{"df":2,"docs":{"1":{"tf":1.0},"3":{"tf":1.0}}}},"df":0,"docs":{}}},"o":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"b":{"df":0,"docs":{},"l":{"df":1,"docs":{"0":{"tf":1.0}}}},"df":0,"docs":{}}}}},"r":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"f":{"df":0,"docs":{},"i":{"df":0,"docs":{},"x":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"i":{"df":0,"docs":{},"m":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"i":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"df":0,"docs":{}}},"o":{"c":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":2,"docs":{"1":{"tf":1.0},"2":{"tf":1.4142135623730951}}}}}},"df":0,"docs":{},"j":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"t":{"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}}},"p":{"df":0,"docs":{},"o":{"df":0,"docs":{},"s":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}}},"u":{"df":0,"docs":{},"l":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":2.0}}}}},"y":{"df":0,"docs":{},"t":{"df":0,"docs":{},"h":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"3":{"df":1,"docs":{"2":{"tf":1.0}}},"df":0,"docs":{}}}}}}},"q":{"df":0,"docs":{},"u":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"3":{"tf":1.4142135623730951}}}}}}},"t":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"3":{"tf":1.0}}}}}}}}},"r":{"df":0,"docs":{},"e":{"a":{"c":{"df":0,"docs":{},"h":{"df":1,"docs":{"1":{"tf":1.0}}}},"d":{"df":0,"docs":{},"i":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"df":0,"docs":{}},"df":0,"docs":{},"f":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}},"l":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"v":{"df":1,"docs":{"1":{"tf":1.0}}}}},"m":{"df":0,"docs":{},"o":{"df":0,"docs":{},"v":{"df":1,"docs":{"1":{"tf":1.0}}}}},"n":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"d":{"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":2,"docs":{"1":{"tf":1.4142135623730951},"2":{"tf":1.4142135623730951}}}}},"df":0,"docs":{}},"p":{"df":0,"docs":{},"l":{"a":{"c":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"o":{"df":2,"docs":{"1":{"tf":1.4142135623730951},"2":{"tf":1.0}},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":0,"docs":{},"o":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.4142135623730951}}}}}}}}}},"q":{"df":0,"docs":{},"u":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"'":{"df":1,"docs":{"1":{"tf":1.0}}},"df":2,"docs":{"0":{"tf":1.7320508075688772},"1":{"tf":1.7320508075688772}}}}}}},"s":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"c":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"s":{"/":{"0":{"0":{"0":{"0":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"<":{"df":0,"docs":{},"f":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{".":{"df":0,"docs":{},"m":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}},"v":{"df":0,"docs":{},"i":{"df":0,"docs":{},"e":{"df":0,"docs":{},"w":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}}},"f":{"c":{"df":3,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":2.6457513110645907},"2":{"tf":1.0}}},"df":0,"docs":{}},"n":{"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"q":{"df":1,"docs":{"1":{"tf":1.7320508075688772}}}}}},"df":0,"docs":{}}},"s":{"c":{"df":0,"docs":{},"o":{"df":0,"docs":{},"p":{"df":0,"docs":{},"e":{"df":1,"docs":{"0":{"tf":1.4142135623730951}}}}}},"df":0,"docs":{},"e":{"df":0,"docs":{},"e":{"df":2,"docs":{"1":{"tf":1.0},"2":{"tf":1.0}}}},"t":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":2,"docs":{"0":{"tf":1.4142135623730951},"1":{"tf":1.4142135623730951}},"e":{"df":0,"docs":{},"p":{"df":1,"docs":{"2":{"tf":1.0}}}},"i":{"df":0,"docs":{},"l":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"u":{"b":{"df":0,"docs":{},"j":{"df":0,"docs":{},"e":{"c":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}},"df":0,"docs":{},"p":{"df":0,"docs":{},"p":{"df":0,"docs":{},"o":{"df":0,"docs":{},"r":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}},"t":{"a":{"df":0,"docs":{},"g":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{},"m":{"df":0,"docs":{},"p":{"df":0,"docs":{},"l":{"a":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}},"e":{".":{"df":0,"docs":{},"m":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}}},"df":0,"docs":{}}}},"x":{"df":0,"docs":{},"t":{"/":{"0":{"0":{"0":{"0":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"1":{"2":{"3":{"4":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{},"r":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"q":{"df":1,"docs":{"1":{"tf":1.0}}}}}},"df":0,"docs":{}}}},"df":1,"docs":{"1":{"tf":1.0}}}}},"i":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.4142135623730951}}}},"t":{"df":0,"docs":{},"l":{"df":1,"docs":{"1":{"tf":1.0}}}}},"o":{"df":0,"docs":{},"w":{"a":{"df":0,"docs":{},"r":{"d":{"df":1,"docs":{"0":{"tf":1.0}}},"df":0,"docs":{}}},"df":0,"docs":{}}},"r":{"a":{"df":0,"docs":{},"n":{"df":0,"docs":{},"s":{"df":0,"docs":{},"p":{"a":{"df":0,"docs":{},"r":{"df":1,"docs":{"0":{"tf":1.0}}}},"df":0,"docs":{}}}}},"df":0,"docs":{}}},"u":{"df":0,"docs":{},"p":{"df":1,"docs":{"1":{"tf":1.0}}},"s":{"df":1,"docs":{"1":{"tf":1.0}},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{">":{"/":{"<":{"df":0,"docs":{},"f":{"df":0,"docs":{},"e":{"a":{"df":0,"docs":{},"t":{"df":0,"docs":{},"u":{"df":0,"docs":{},"r":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"a":{"df":0,"docs":{},"m":{"df":1,"docs":{"1":{"tf":1.0}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{}}}},"df":0,"docs":{}},"df":0,"docs":{}},"df":0,"docs":{}}}},"df":0,"docs":{}}}}}},"v":{"2":{".":{"0":{".":{"df":0,"docs":{},"m":{"d":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}},"df":0,"docs":{}}},"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{},"e":{"df":0,"docs":{},"r":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"2":{"tf":1.0}}}}}}}},"i":{"df":0,"docs":{},"s":{"df":0,"docs":{},"i":{"df":0,"docs":{},"t":{"df":1,"docs":{"2":{"tf":1.0}}}}}}},"w":{"a":{"df":0,"docs":{},"r":{"df":0,"docs":{},"r":{"a":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"0":{"tf":1.0}}}}},"df":0,"docs":{}}}},"df":0,"docs":{},"i":{"df":0,"docs":{},"p":{"df":1,"docs":{"1":{"tf":1.4142135623730951}}}},"o":{"df":0,"docs":{},"r":{"df":0,"docs":{},"k":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":1.0}},"f":{"df":0,"docs":{},"l":{"df":0,"docs":{},"o":{"df":0,"docs":{},"w":{"df":2,"docs":{"0":{"tf":1.0},"1":{"tf":2.0}}}}}}}}}},"y":{"df":0,"docs":{},"o":{"df":0,"docs":{},"u":{"'":{"d":{"df":1,"docs":{"1":{"tf":1.0}}},"df":0,"docs":{}},"df":0,"docs":{},"r":{"df":0,"docs":{},"s":{"df":0,"docs":{},"e":{"df":0,"docs":{},"l":{"df":0,"docs":{},"f":{"df":1,"docs":{"1":{"tf":1.0}}}}}}}}}}}},"title":{"root":{"b":{"df":0,"docs":{},"u":{"df":0,"docs":{},"i":{"df":0,"docs":{},"l":{"d":{"df":1,"docs":{"2":{"tf":1.0}}},"df":0,"docs":{}}}}},"c":{"df":0,"docs":{},"o":{"df":0,"docs":{},"m":{"df":0,"docs":{},"m":{"df":0,"docs":{},"e":{"df":0,"docs":{},"n":{"df":0,"docs":{},"t":{"df":1,"docs":{"0":{"tf":1.0}}}}}}},"n":{"df":0,"docs":{},"t":{"df":0,"docs":{},"r":{"df":0,"docs":{},"i":{"b":{"df":0,"docs":{},"u":{"df":0,"docs":{},"t":{"df":1,"docs":{"1":{"tf":1.0}}}}},"df":0,"docs":{}}}}}}},"df":0,"docs":{},"p":{"df":0,"docs":{},"r":{"df":0,"docs":{},"o":{"c":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"s":{"df":1,"docs":{"2":{"tf":1.0}}}}}},"df":0,"docs":{}}}},"q":{"df":0,"docs":{},"u":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":0,"docs":{},"i":{"df":0,"docs":{},"o":{"df":0,"docs":{},"n":{"df":1,"docs":{"3":{"tf":1.0}}}}}}}}}},"r":{"df":0,"docs":{},"e":{"df":0,"docs":{},"q":{"df":0,"docs":{},"u":{"df":0,"docs":{},"e":{"df":0,"docs":{},"s":{"df":0,"docs":{},"t":{"df":1,"docs":{"0":{"tf":1.0}}}}}}}}}}}},"pipeline":["trimmer","stopWordFilter","stemmer"],"ref":"id","version":"0.9.5"},"results_options":{"limit_results":30,"teaser_word_count":30},"search_options":{"bool":"OR","expand":true,"fields":{"body":{"boost":1},"breadcrumbs":{"boost":1},"title":{"boost":2}}}});