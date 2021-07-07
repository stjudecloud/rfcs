#!/usr/bin/env bash
set -e

function build() {
  _BRANCH=$1;
  _DIR=$2;
  _DRAFTS_FILE=$3;

  echo "  [*] Building ${_BRANCH} to ${_DIR}"
  git checkout "${_BRANCH}"
  
  [[ -d src/ ]] && rm -rf src/
  mkdir -p ${_DIR} src/

  printf "# Summary\n\n[Introduction](introduction.md)\n\n" > src/SUMMARY.md
  cp README.md src/introduction.md

  for RFC_FILE in $(ls text/* | sort)
  do
      echo "- [$(basename ${RFC_FILE} ".md")]($(basename ${RFC_FILE}))" >> src/SUMMARY.md
      cp ${RFC_FILE} src
  done

  printf "\n[RFC Drafts]($(basename ${_DRAFTS_FILE}))\n" >> src/SUMMARY.md
  cp ${_DRAFTS_FILE} src
  mdbook build -d "${BRANCH_DIR}"
  cp -R resources "${BRANCH_DIR}"
  rm -rf src/
}

# Stash starting branch name
STARTING_BRANCH_NAME="$(git rev-parse --abbrev-ref HEAD)"

# Create book directory
BOOK_DIR="./book"
mkdir -p "$BOOK_DIR"

DRAFTS_FILE="drafts.md"

# Get branch list
git fetch --all
BRANCHES=()
while IFS= read -r line; do
    BRANCHES+=( "$line" )
done < <( git branch --list --all | sed 's,\*,,g' | xargs -n1 | grep "remotes/origin" | sed 's,remotes/origin/,,g' | sort | uniq | grep -e 'master' -e 'rfcs/' )

echo "== Creating Drafts File =="
printf "# Drafts\n\n" > "${DRAFTS_FILE}"
echo "The following are _candidate_ RFCs that are being rendered for easy review. They are *not* accepted St. Jude Cloud RFCs. For more information please see [the associated pull request](https://github.com/stjudecloud/rfcs/pulls)." >> "${DRAFTS_FILE}"
printf "\n\n" >> "${DRAFTS_FILE}"

for CURRENT_BRANCH in "${BRANCHES[@]}"; do
  if [[ "${CURRENT_BRANCH}" != "master" ]]; then
    echo "- [${CURRENT_BRANCH}](https://stjudecloud.github.io/rfcs/branches/$CURRENT_BRANCH)" >> "${DRAFTS_FILE}"
  fi
done

# Loop through rfc branches, build and copy each to output dir
echo "== Build RFC Branches =="
for CURRENT_BRANCH in "${BRANCHES[@]}"; do 
  BRANCH_DIR="${BOOK_DIR}/branches/${CURRENT_BRANCH}"
  if [[ "${CURRENT_BRANCH}" == "master" ]]; then
    BRANCH_DIR="${BOOK_DIR}/"
  fi
  build "${CURRENT_BRANCH}" "${BRANCH_DIR}" "${DRAFTS_FILE}"
done

rm ${DRAFTS_FILE}
git checkout "${STARTING_BRANCH_NAME}"
