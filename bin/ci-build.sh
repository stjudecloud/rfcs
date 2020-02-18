#!/usr/bin/env bash
set -e

# Get the branch name
BRANCH_NAME="$(git rev-parse --abbrev-ref HEAD)"
FINAL_DIRECTORY="book/branches/${BRANCH_NAME}/"
if [[ "${BRANCH_NAME}" == "master" ]]; then
  FINAL_DIRECTORY="book/"
fi

mkdir -p $FINAL_DIRECTORY

[[ ! -d src ]] && mkdir src
printf "# Summary\n\n[Introduction](introduction.md)\n" > src/SUMMARY.md

for f in $(ls text/* | sort)
do
    echo "- [$(basename $f ".md")]($(basename $f))" >> src/SUMMARY.md
    cp $f src
done

cp README.md src/introduction.md

mdbook build -d "$FINAL_DIRECTORY"
cp -R resources "$FINAL_DIRECTORY"