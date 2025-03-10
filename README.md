# fur_hotspot_mutations

This repository contians code and analysis results for identifying artefacts and rescuing false negative variant calls in the FUR Felis Catus dataset.

|                         Main                         |                         Develop                          |
| :----------------------------------------------------: | :------------------------------------------------------: |
| [![pipeline status][main-pipe-badge]][main-branch] | [![pipeline status][develop-pipe-badge]][develop-branch] |

[main-pipe-badge]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/badges/main/pipeline.svg
[main-branch]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/-/commits/main
[develop-pipe-badge]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/badges/develop/pipeline.svg
[develop-branch]: https://gitlab.internal.sanger.ac.uk/DERMATLAS/fur/fur_hotspot_mutations/-/commits/develop


## Table of Contents
- [Installation and Usage](#installation-and-usage)
- [Directory Structure](#directory-structure)
- [CI/CD Pipeline Overview](#cicd-pipeline-overview)
- [Development Installation](#development-installation)
  - [Local/Remote Machine Setup](#localremote-machine-setup)
    - [Prerequisites](#prerequisites)
    - [Setup Steps](#setup-steps)
  - [VSCode Devcontainer Setup](#vscode-devcontainer-setup)
    - [Prerequisites](#prerequisites-1)
    - [Setup Steps](#setup-steps-1)
  - [Building with Docker](#building-with-docker)
- [Running Tests](#running-tests)
- [Rolling Out a Tagged Release](#rolling-out-a-tagged-release)
- [Publishing package to PyPi](#publishing-package-to-pypi)

## Installation and Usage

This package is hosted on the [internal Gitlab PyPi server](https://gitlab.internal.sanger.ac.uk/groups/team113sanger/-/packages/). To install it, run the following command:

```bash
# If using poetry
PROJECT_URL="https://gitlab.internal.sanger.ac.uk/api/v4/projects/5569/packages/pypi/simple"
poetry source add --priority=supplemental gitlab-pypi $PROJECT_URL
poetry add fur_hotspot_mutations --source gitlab-pypi

# If using pip
PROJECT_URL="https://gitlab.internal.sanger.ac.uk/api/v4/projects/5569/packages/pypi/simple"
pip install fur_hotspot_mutations --index-url $PROJECT_URL

# If using poetry but pulling the pacakge from the group level, you will need a personal access token or deploy token
GROUP_URL="https://gitlab.internal.sanger.ac.uk/api/v4/groups/882/-/packages/pypi/simple"
poetry source add --priority=supplemental gitlab-pypi-grp $GROUP_URL
poetry config http-basic.gitlab-pypi-grp <your_username> <your_token>
poetry add fur_hotspot_mutations --source gitlab-pypi-grp
```

If you see **authentication errors with** `pip` or `poetry`, your Gitlab repo has set it's PyPI to private - you can either use Personal Access Tokens or Deploy Tokens to authenticate. Or more conveniently, you can **change the project's PyPI settings to public**. Visit the Gitlab Repo's settings and you will want `Settings > General > Visibility, project features, permissions > Package Registry > Allow anyone to pull from Package Registry` to be set to `Enabled` and then **save the changes**.

Once the `fur_hotspot_mutations` package is installed, it provides two command-line scripts: `add` and `multiply`. These scripts can be called directly from the command line.

## Directory structure

```
.
├── .devcontainer
│   └── devcontainer.json    # Used by VSCode to create a development container
├── cicd_scripts/
│   ├── build_package_and_publish.sh       # Publishes the package to the Gitlab PyPi registry
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

The pipeline makes use of predefined CI/CD templates included from the [team113sanger/common/cicd-template](https://gitlab.internal.sanger.ac.uk/team113sanger/common/cicd-template) repository. These components provide predefined variables, scripts, and rules used in the pipeline. They are version controlled and included into the `.gitlab-ci.yml` file using the `include` keyword. For more details on each component, please refer to the repo.

### Stages

#### Build

The `build` job is responsible for building the Docker image for the project. It uses Docker's BuildKit to speed up the build process. This job is run on every commit, excluding those tagged for release.

It pushes a **candidate Docker image** for user in subsequent jobs as it is not trusted yet.

#### Test

The `test` stage includes the `unit-test` and `format-test` job, which runs unit tests on the project using the **candidate Docker image** built in the `build` stage. This job is also run on every commit, excluding those tagged for release.

#### Publish

The `publish` stage includes two jobs:

- `publish-develop`: This job publishes the **candidate Docker image** to the gitlab Docker registry, tagging it `develop-branch-latest`. It only runs on commits to the `develop` branch.
- `publish-tagged_and_latest`: This job publishes the **candidate Docker image**, tagging it `latest` and `$SEM_VER` i.e. the git tag version. It only runs on tagged commits.

## Development Installation

This section provides a guide on how to set up your development environment to
contribute to this project. We use Poetry for dependency management
and provide two different setup options: *Local/Remote setup* and *VSCode
devcontainer setup*.

### Local/Remote Machine Setup

#### Prerequisites:
 -  Python 3.9 (or later) installed
    - (if on the farm, do `module load python/3.9.2` or similar)
 - Poetry installed
    - if on the farm, do `module load /software/team113/modules/modulefiles/poetry/1.8.2`
    - if on the farm and can't access the module, create a virtual env and install `poetry` (will be covered below)
    - if on you are on OpenStack or local machine, follow the [official instructions here](https://python-poetry.org/docs/#installation)
- Git + hubflow installed
    - if on the farm, do `module load /software/CASM/modules/modulefiles/git/2.21.0`
    - if on your own machine or OpenStack, follow the [official instructions here](https://datasift.github.io/gitflow/TheHubFlowTools.html)

#### Setup steps:
1. `git clone` this repo and `cd` into it
1. Create a python virtual environment (venv) and activate it
    - `python -m venv .venv`
    - `source .venv/bin/activate`
1. Optional, install poetry in the venv if you can't access the module on the farm
    - `pip install poetry==1.8.2`
1. Once you are in a virtual env, do `poetry install` to install all the Python dependencies.
   * If you experience any of the following exceptions `PromptDismissedException`, `ItemNotFoundException`, `DBusErrorResponse`, `InitError`, there may be an issue with the configuration of keyrings in Linux. This is a Poetry ^1.4.0 issue. To fix it:
     * `export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring` then `poetry install`
1. Do `git hf init` to intialise the local version of the repo with *Hubflow*
1. Do `pre-commit install --install-hooks` to install your local git hooks (linting, formatting)

### VSCode Devcontainer Setup

This setup allows you to develop this project inside a Docker container using
Visual Studio Code's Dev Container feature. This setup provides a **consistent** and
fully featured development environment, which includes all the necessary tools
and packages.

**It cannot be understated that this setup is very consistent but cannot be done
on the farm**. It is recommended to use this setup on your own machine or on
OpenStack.

The files which enable you to work from within an container can be found in
`.devcontainer/devcontainer.json`. This turns the Docker image into a
development Docker image.

#### Prerequisites:
- Docker installed (either locally or on OpenStack)
- Visual Studio Code installed
  - You can download it [here](https://code.visualstudio.com/download).
- "Remote - Containers" extension installed in VSCode
  - You can install the extension directly from the [Visual Studio Marketplace](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) or by searching for "Remote Containers" in the VSCode Extensions view

#### Setup steps:
1. `git clone` this repo and `cd` into it
1. Open the project in VSCode.
1. Start the Dev Container.

   Click on the green "open a remote window" button in the lower left corner of VSCode, or use the `F1` command palette and select "Remote-Containers: Reopen in Container". This will start building the Docker container as described in the project's `.devcontainer/devcontainer.json` file. The build process can take a few minutes when run for the first time.

1. Wait for the Dev Container to be built.

   VSCode will use the configuration in the `.devcontainer` folder of the project to build a Docker container, install all necessary software, and open a VSCode window connected to that container.

1. Setup `git` in the container.

   The container is a fresh environment, so you'll need to configure `git` to use your name and email address. You can do this by running the following commands in the terminal in VSCode:

   ```bash
   git config --global user.name "Your Name"
   git config --global user.email "your@email.com"

   # You can check git is work in the container by running
   git fetch
   ```

   **IMPORTANT**: The SSH keys you use to authenticate with Gitlab are proxied into the container through a Linux socket, so you don't need to do anything else to authenticate with Gitlab... but [troubleshooting advice is here](https://code.visualstudio.com/remote/advancedcontainers/sharing-git-credentials). Where the instructions mention `github_rsa` it may just be `id_rsa` for you.

   **IMPORTANT**: The ssh public key that is sent through the socket comes from your **LOCAL MACHINE** i.e. your Sanger laptop and not the remote OpenStack host. This is not obvious. You will need to add your Sanger laptop's public key to your Gitlab account. You can find your public key by running `cat ~/.ssh/id_rsa.pub` on your Sanger laptop. It is also important to ensure the ssh key-agent is running on your Sanger laptop. You can check this by running `eval "$(ssh-agent -s)"` on your Sanger laptop; otherwise you will need to update you ssh config to start the agent on startup or your bashrc.

1. Start Developing!

   Once VSCode connects to the container, you're in the development environment defined by the container! The terminal in VSCode now runs inside the container, and you should have all the tools you need for development pre-installed: `hubflow`, `poetry`, `pre-commit`, and all the Python dependencies. You can even commit from within the container!

Please check the [official VSCode documentation](https://code.visualstudio.com/docs/devcontainers/containers) for more detailed information on developing inside a container with VSCode.

### Building with Docker

While VSCode is the recommended way to develop this project as it handle the setup of the Docker container and the installation of all the dependencies, you can also build the Docker container manually.

```bash
# Build the Docker image (buildkit is required), without using the cache and targeting the base_stage only
export DOCKER_BUILDKIT=1
MY_IMAGE="fur_hotspot_mutations:local"
docker build -t $MY_IMAGE --target base_stage --no-cache .

# Run the Docker container
docker run --rm $MY_IMAGE add 4.5 5.6 # example command
docker run -it --rm $MY_IMAGE bash # to enter an interactive bash shell in the container
```

If you use docker-compose, you can also build and run the container indefinitely, with the src & test directories mounted as volumes:

```bash
# Build & run the Docker image in the background
docker compose up --build -d
docker container ls # to find the container name
docker exec -it <container_name> bash # to enter an interactive bash shell in the container
```


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
