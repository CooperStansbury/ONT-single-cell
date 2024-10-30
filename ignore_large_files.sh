#!/bin/bash

# --- Parameters ---
DIRECTORY_TO_SEARCH="${1:-.}"  # Defaults to the current directory
GITIGNORE_PATH="${2:-.gitignore}"  # Defaults to .gitignore in the current directory
MAX_FILE_SIZE_MB="${3:-100}"  # Defaults to 100 MB

# --- Print script parameters ---
echo "Searching in directory: ${DIRECTORY_TO_SEARCH}"
echo "Using .gitignore file: ${GITIGNORE_PATH}"
echo "Maximum file size (MB): ${MAX_FILE_SIZE_MB}"

# --- Find files exceeding the size limit and add them to .gitignore ---
find "${DIRECTORY_TO_SEARCH}" -type f -size +"$((MAX_FILE_SIZE_MB * 1024 * 1024))"c -print0 | while IFS= read -r -d '' FILE; do
    echo "Adding ${FILE} to ${GITIGNORE_PATH}"
    echo "${FILE}" >> "${GITIGNORE_PATH}"
done