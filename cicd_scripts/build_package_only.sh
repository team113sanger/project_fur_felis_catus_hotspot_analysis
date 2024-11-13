#!/bin/bash

set -e

THIS_PROGRAM=./cicd_scripts/build_package_only.sh

# Check if pyproject.toml exists in the current working directory
if [[ ! -f "./pyproject.toml" ]]; then
  >&2 echo "Error: Script is designed to be run while CD to a directory containing pyproject.toml"
  >&2 echo "Error: Hence the typical invocation is 'cd <project_dir> && ${THIS_PROGRAM}'"
  exit 1
fi

# Build using poetry
BUILD_COMMAND="poetry build --no-ansi --format wheel"
echo "Running command: ${BUILD_COMMAND}" >&2
POETRY_BUILD_OUTPUT=$(eval "${BUILD_COMMAND}")

# Echo the output of the build command to stderr
echo "$POETRY_BUILD_OUTPUT" >&2

# Parse the input string to get the wheel file name
# Awk will split the string by spaces and print the third element
# For example: '- Built example_python_cicd-0.1.0-py3-none-any.whl' will be split into
# '-', 'Built', 'example_python_cicd-0.1.0-py3-none-any.whl'
WHEEL_FILE=$(echo "$POETRY_BUILD_OUTPUT" | grep 'Built' | grep '.whl' | awk '{print $3}')

# Define the path to the dist directory
DIST_DIR="dist"

# Check if the wheel file exists in the dist directory
if [[ -e "${DIST_DIR}/${WHEEL_FILE}" ]]; then
  echo "${DIST_DIR}/${WHEEL_FILE}"
else
  # If the file doesn't exist, write an error message to stderr and exit with status code 1
  >&2 echo "Error: ${WHEEL_FILE} does not exist in ${DIST_DIR}"
  exit 1
fi
