#!/usr/bin/env bash
set -e

[[ ! -d src ]] && mkdir src
echo "[Introduction](introduction.md)\n" > src/SUMMARY.md

for f in $(ls text/* | sort)
do
    echo "- [$(basename $f ".md")]($(basename $f))" >> src/SUMMARY.md
    cp $f src
done

cp README.md src/introduction.md
mdbook build
