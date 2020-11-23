#!/usr/bin/env bash

[[ -d book/ ]] && rm -rf book/
[[ -d src/ ]] && rm -rf src/
mkdir -p src/

printf "# Summary\n\n[Introduction](introduction.md)\n\n" >src/SUMMARY.md
cp README.md src/introduction.md

for RFC_FILE in $(ls text/* | sort); do
    echo "- [$(basename ${RFC_FILE} ".md")]($(basename ${RFC_FILE}))" >>src/SUMMARY.md
    cp ${RFC_FILE} src
done

mdbook build
