#!/bin/bash

set -e

# Check for necessary environment variables
if [[ -z "${GITLAB_DEPLOY_USERNAME}" ]]; then
  echo "Error: GITLAB_DEPLOY_USERNAME is not defined. Please define GITLAB_DEPLOY_USERNAME in the environment variables." >&2
  echo "Deploy tokens can be created in the GitLab group's settings > repository > deploy token." >&2
  echo "If using the GitLab CICD to run this script you could use value 'gitlab-ci-token' as the deploy username." >&2
  exit 1
fi

if [[ -z "${GITLAB_DEPLOY_TOKEN}" ]]; then
  echo "Error: GITLAB_DEPLOY_TOKEN is not defined. Please define GITLAB_DEPLOY_TOKEN in the environment variables." >&2
  echo "Deploy tokens can be created in the GitLab group's settings > repository > deploy token." >&2
  echo "If using the GitLab CICD to run this script you could use the value '\$CI_JOB_TOKEN' for the deploy token." >&2
  exit 1
fi

# Define the repository URL and name
# Url is the Sanger internal GitLab API url
PROJECT_ID="5569" # This is the project ID found in the GitLab repository settings
REPOSITORY_NAME="gitlab-sanger-$PROJECT_ID" # This name is arbitrary
REPO_URL="https://gitlab.internal.sanger.ac.uk/api/v4/projects/${PROJECT_ID}/packages/pypi/"

# Register the repository using poetry
echo "Poetry registering repository ${REPOSITORY_NAME} with url ${REPO_URL}"
poetry config --unset "repositories.$REPOSITORY_NAME" &>/dev/null || true
poetry config "repositories.$REPOSITORY_NAME" "${REPO_URL}"

# Clearing the dist directory
echo "Clearing the dist directory, to allow Poetry to build the package"
rm -rf dist

echo "Building package using poetry"
poetry build

# Publish the package using poetry
# --username and --password are used for authentication
# If the package already exists on the repository, the publish command will fail
# but we want tolerate this error without causing the script to error out
PUBLISH_COMMAND="poetry publish  --repository $REPOSITORY_NAME --username ${GITLAB_DEPLOY_USERNAME} --password ${GITLAB_DEPLOY_TOKEN} --no-interaction"
echo "Running publish command: ${PUBLISH_COMMAND}"


set +e  # Temporarily disable exit on error, to collect the output and exit code
# Execute the command, and then immediately capture the status and stderr
TEMPFILE_STDERR=$(mktemp)
eval "${PUBLISH_COMMAND}" 2> $TEMPFILE_STDERR
PUBLISH_COMMAND_EXIT_CODE=$?
PUBLISH_COMMAND_STDERR=$(cat $TEMPFILE_STDERR)
rm ${TEMPFILE_STDERR}
set -e  # Re-enable exit on error


if [ $PUBLISH_COMMAND_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "ERROR: Command: ${PUBLISH_COMMAND}"
    echo "ERROR: Exit code: ${PUBLISH_COMMAND_EXIT_CODE}"
    echo "ERROR: Stderr:"
    echo "$PUBLISH_COMMAND_STDERR"
    if [[ $PUBLISH_COMMAND_STDERR == *"File name has already been taken"* ]]; then
        echo "ERROR: Known 400 HTTP error - file name has already been taken."
        echo "INFO: Package already exists on Gitlab PyPi, so cannot be re-published (error handled)."
    elif [[ $PUBLISH_COMMAND_STDERR == *"Description is too long"* ]]; then
        echo "ERROR: Known 400 HTTP error - GitLab server currently cannot process pyproject.toml with a large README."
        echo "INFO: See https://gitlab.com/gitlab-org/gitlab/-/issues/431505 for more info."
        echo "INFO: Solution? Update GitLab to 16.6 or newer, or remove the 'readme' key from pyproject.toml and try again."
    else
        echo "ERROR: Unknown error when publishing package to Gitlab PyPi."
    fi
    exit 1
else
    echo "Package published successfully to Gitlab PyPi."
fi

echo "Done"
