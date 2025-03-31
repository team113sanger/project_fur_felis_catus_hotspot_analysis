# fur_hotspot_mutations

|                         Main                         |                         Develop                          |
| :----------------------------------------------------: | :------------------------------------------------------: |
| [![pipeline status][master-pipe-badge]][master-branch] | [![pipeline status][develop-pipe-badge]][develop-branch] |

[master-pipe-badge]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/badges/master/pipeline.svg
[master-branch]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/-/commits/master
[develop-pipe-badge]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/badges/develop/pipeline.svg
[develop-branch]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/-/commits/develop

## Usage

Within the Docker image (see below) you can run the following command:

```
usage: fur_hotspot_mutations [-h] [--version] COMMAND ...

A toolkit for identifying hotspot mutations in tumor cohorts and rescuing false-negative variant calls from Caveman and Pindel.
The program processes MAF files to evaluate mutations against specific criteria, rescue false-negatives, remove false-positives, and reclassify suspected germline mutations.

positional arguments:
  COMMAND
    extract_hotspot_mutations
                        Identify and extract hotspot mutations from a MAF file.
    mpileup_variant_filter
                        Filter variant calls using mpileup data.
    update_maf_variants
                        Add or remove variants in a MAF file.
    check_mpileup_samples_from_maf
                        Process MAF and BAM files using samtools mpileup.

options:
  -h, --help            show this help message and exit
```

Please see tests for examples

## Docker Image

This project hosts Docker images on Quay.io. Please see [https://quay.io/repository/team113sanger/fur_hotspot_mutations](https://quay.io/repository/team113sanger/fur_hotspot_mutations?tab=tags).

## Summary

This repository contains code for identifying artifacts and rescuing false-negative variant calls from Caveman and Pindel calls. It also includes analysis results produced from the FUR Felis Catus dataset.

In brief, the tool takes an input MAF file (generated with the Dermatlas processing pipeline) that describes a cohort of tumors and identifies hotspot mutations within the cohort. The tool then generates a read pileup for each hotspot locus in each sample and counts the number of reads attributable to alternative alleles in the pileup.

Hotspot mutations in the MAF are then evaluated against several criteria:

- Is the variant allele frequency (VAF) > 1%?
- Is the number of ALT allele reads in the tumor sample > 5?
- Is the number of ALT allele reads in the normal sample < 3?

When these criteria are met, the mutation is considered a false-negative and new occurrences of the hotspot mutation are added to the MAF file. Conversely, in samples where an identified hotspot mutation fails to meet these criteria, it is considered a false positive variant and removed from the MAF.

The tool also looks at calls at the cohort level. When the number of tumor-normal pairs flagged as false negatives for a variant is greater than 3, the variant is reclassified as a germline mutation and removed.


## Table of Contents
- [Installation and Usage](#installation-and-usage)
- [Directory Structure](#directory-structure)

## Directory structure

```
.
├── .devcontainer
│   └── devcontainer.json    # Used by VSCode to create a development container
├── cicd_scripts/
│   ├── build_package_and_publish.sh       # Publishes the package to the Sanger Gitlab PyPi registry
│   └── format_test.sh       # CI/CD script to run formatting and linting
├── .dockerignore
├── .gitignore
├── .gitlab-ci.yml           # Gitlab CI/CD config file
├── .pre-commit-config.yaml  # Pre-commit config file
├── Dockerfile
├── README.md
├── docker-compose.yml
├── pyproject.toml           # Python package config file (dependencies add here)
├── poetry.lock              # Poetry dependency lock file (generated)
├── poetry.toml              # Poetry config file
├── requirements.txt         # Alternative package dependency file (generated)
├── src
│   └── fur_hotspot_mutations  # The python package of this repo (src code)
└── tests
    ├── integration          # Integration tests go here
    └── unit                 # Unit tests go here
```

## CI/CD Pipeline Overview

The project's CI/CD pipeline is configured in `.gitlab-ci.yml` and comprises three stages: `build`, `test`, and `publish`.

### Local/Remote Machine Setup

#### Prerequisites:
 -  Python 3.9 (or later) installed
 - Poetry installed
    - Follow the [official instructions here](https://python-poetry.org/docs/#installation)
- Git + hubflow installed

## Running Tests

This project uses pytest as the main testing framework and all tests are located in the `tests` directory. For formatting and linting, we use `pre-commit` hooks which are installed in the setup steps above.

### Unit & Integration Tests

Before running the tests, ensure that your virtual environment is activated and all the required packages are installed. If you followed the above setup steps, your environment should be ready.

To run the tests, navigate to the project's root directory and execute the following command:

```bash
# To run all tests
pytest
# If you want to see a more verbose output
pytest -sv
# If you want to have the tests stop after the second failure
pytest --maxfail 2
# If you want to run only the failing tests from the last run
pytest --lf
```

### Pre-commit (Linting & Formatting)

This project uses `pre-commit` hooks to ensure that all code is formatted and linted before committing. That means whenever you commit the code will be checked and then check again on the remote GitLab server. If the code is not formatted or linted, the commit will fail.

To run the hooks manually, execute the following command:

```bash
pre-commit run --all-files
# Or use the CICD compatible script
./cicd_scripts/format_test.sh
```

## Rolling out a Tagged Release
This repo the GitFlow branching model and uses [hubflow](https://datasift.github.io/gitflow/TheHubFlowTools.html) as a tool to enable this from the CLI.

```bash
# Switch to the develop branch
git checkout develop

# Start a new release branch e.g. 0.1.0 not v0.1.0
git hf release start <project_version>
```

Now, do the following things:
* `CHANGELOG.md`: Under the heading of the newest release version, describe what was changed, fixed, added.
* `pyproject.toml`: Increment the project version to the current release version
* Commit these changes

Finally

```bash
git hf release finish <project_version>
```

### Publishing package to PyPi
This repo uses [poetry](https://python-poetry.org/) to manage dependencies and packaging. If pushing to the intenral GitLab pypi, you can use a [group token](https://gitlab.internal.sanger.ac.uk/help/user/packages/package_registry/index#authenticate-with-the-registry) instead of a personal access token.

To publish a new version to PyPi, do the following:

```bash
export GITLAB_DEPLOY_USERNAME=<your_pypi_token>
export GITLAB_DEPLOY_TOKEN=<your_pypi_password>
./cicd_scripts/build_package_and_publish.sh
```

To build a package without publishing, do the following:

```bash
poetry build
# Or if you need a script that returns the path to a Python wheel file
./cicd_scripts/build_package_only.sh
```
