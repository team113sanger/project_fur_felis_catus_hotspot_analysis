# syntax=docker/dockerfile:1
# MUTLI-STAGE DOCKERFILE for building a Python image
# STAGE 1: base_stage is the base image for both the development and dev-staging/staging/production images
# STAGE 2: development_only_stage is the image used for development (optimised for VSCode)

########################
# STAGE 1: base_stage #
########################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM python:3.11-slim-bullseye as base_stage

# Set environment variables.
# 1. Force Python stdout and stderr streams to be unbuffered.
# 2. Do not build if pip can't connect to the internet
# 3. Do not cache pip installs
# 4. Set the poetry version which will handle the rest of the depenedency installation
# 5. The mode in which Django will be executed
USER root
ENV \
    PYTHONUNBUFFERED=1 \
    PIP_DEFAULT_TIMEOUT=100 \
    PIP_NO_CACHE_DIR=1 \
    PIP_VERSION=23.3.1 \
    POETRY_VERSION=1.8.2 \
    DATA_DIRECTORY="/data" \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    POETRY_VIRTUALENVS_CREATE=false \
    POETRY_VIRTUALENVS_IN_PROJECT=false

ENV \
    USER_BASHRC="${USER_DIRECTORY}/.bashrc" \
    USER_BIN_DIRECTORY="${USER_DIRECTORY}/.local/bin" \
    SSH_DIR="${USER_DIRECTORY}/.ssh" \
    VENV_DIRECTORY="${OPT_DIRECTORY}/venv" \
    POETRY_HOME="${OPT_DIRECTORY}/poetry" \
    PIPX_HOME="${OPT_DIRECTORY}/pipx" \
    PIPX_BIN_DIR="${OPT_DIRECTORY}/pipx/bin" \
    POETRY_CACHE_DIR="${OPT_DIRECTORY}/poetry-cache" \
    PROJECT_DIRECTORY="${OPT_DIRECTORY}/repo" \
    LOGGING_DIRECTORY="${DATA_DIRECTORY}/logs"

ENV PATH=${POETRY_HOME}/bin:${PIPX_BIN_DIR}:${PATH}

RUN \
    useradd "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && mkdir -p "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}" "${OPT_DIRECTORY}" "${POETRY_CACHE_DIR}" "${PIPX_HOME}" "${POETRY_HOME}" "${VENV_DIRECTORY}" "${PIPX_BIN_DIR}"\
    && chown -R "${USER_NAME}:${USER_NAME}" "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}" "${USER_DIRECTORY}" "${OPT_DIRECTORY}" \
    && chmod -R 755 "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}" "${USER_DIRECTORY}" "${OPT_DIRECTORY}"


# Update System
# Install system packages required by Python packages
# We include:
# - curl + wget, for debugging or downloading files
# - nano, for debugging & a minimal editor
# - git, so that pre-commit can be installed and run
# - tree, so that the project structure can be viewed in the container
# - netcat, to do network debugging
# - openssh-client, to be able to do git operations over ssh & also scp for backup-and-restore operations
# - build-essential, meta-packages that are essential to compile software including gcc, g++, make, etc.
# - pkg-config, to manage compile and link flags for libraries
# BuildKit logic to cache apt packages
# - https://vsupalov.com/buildkit-cache-mount-dockerfile/
# - https://github.com/moby/buildkit/blob/master/frontend/dockerfile/docs/reference.md#run---mounttypecache
# - https://stackoverflow.com/a/72851168
RUN rm -f /etc/apt/apt.conf.d/docker-clean; echo 'Binary::apt::APT::Keep-Downloaded-Packages "true";' > /etc/apt/apt.conf.d/keep-cache
RUN \
    python3 -m pip install --upgrade pip==${PIP_VERSION} && \
    apt-get update --quiet && \
    apt-get install --yes --quiet --no-install-recommends \
    build-essential \
    pkg-config \
    nano \
    curl \
    wget \
    git \
    tree \
    openssh-client \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m pip install pipx \
    && pipx ensurepath \
    && pipx install poetry==$POETRY_VERSION \
    && pipx inject -f poetry poetry-plugin-export


# As the non-root user we install pipx and poetry so that they are available in
# the command line for subsequent package installations.
USER $USER_NAME


# VENV STEP
# There are two motivations:
# 1. To keep the image small and well configured, we want to re-use system
#    installed packages. This is especially important with tensorflow which is a
#    large package.
# 2. We want the extra packages to be installed in a virtual environment so that
#    when working in a Singularity container which is a read-only file system,
#    we can still install packages (by binding and reinstalling into venv directory)

WORKDIR $OPT_DIRECTORY
RUN \
    python3 -m venv ${VENV_DIRECTORY} &&  \
    echo "source ${VENV_DIRECTORY}/bin/activate" >> "${USER_BASHRC}"
ENV \
    PATH="${VENV_DIRECTORY}/bin:${PATH}" \
    VIRTUAL_ENV=${VENV_DIRECTORY}


# Copy dependency control files and install dependencies
# 1. only the dependency control files for Poetry (equivalent to "requirements.txt"
#    if you have never seen Poetry *.lock and *.toml files before)
#    The unsual square brakets in `poetry.loc[k]` provides GLOB/REGEX syntax which ensures
#    IDEMPOTENCY for the build (i.e. Docker builds correctly the first and every subsequent time).
#    Initially, `poetry.lock` may be absent as it is only created after the first `poetry install`
#    This syntax lets COPY conditionally include it without the build failing if missing.
# 2. .gitignore as it is used by precommit-check in a post build stage of the CICD
WORKDIR $PROJECT_DIRECTORY
COPY --chown="${USER_NAME}:${USER_NAME}" [".gitignore", "pyproject.toml", "poetry.loc[k]", "./"]
RUN \
    poetry install --no-root --no-directory  \
    && chown -R "${USER_NAME}:${USER_NAME}" ${VENV_DIRECTORY}

# Copy the src, test and other code of the project into the container.
# Then install the project itself via Poetry
COPY --chown="${USER_NAME}:${USER_NAME}" . "${PROJECT_DIRECTORY}"
RUN \
    poetry install \
    && chown -R "${USER_NAME}:${USER_NAME}" "${USER_DIRECTORY}" "${PROJECT_DIRECTORY}" "${DATA_DIRECTORY}"

# Use user "admin" to run the build commands below and the server itself.
USER "${USER_NAME}"

###################################
# STAGE 2: development_only_stage #
# - This stage is optional        #
# - It is optimised for VSCode    #
###################################

# To develop from the container we need add some extra directories to work nicely with VSCode
FROM base_stage as development_only_stage
USER root

# Install hubflow to allow for git flow style development & conditional install
# sudo, giving the user passwordless sudo privileges
WORKDIR "${USER_DIRECTORY}"
ARG HAS_SUDO="${HAS_SUDO:-0}"
RUN git config --global --add safe.directory "${USER_DIRECTORY}/gitflow" \
    && git config --global --add safe.directory "${USER_DIRECTORY}/gitflow/shFlags" \
    && git clone https://github.com/datasift/gitflow \
    && chown -R "${USER_NAME}:${USER_NAME}" gitflow \
    && cd gitflow \
    && ./install.sh \
    && if [ "${HAS_SUDO}" = "1" ]; then \
    apt-get update -y \
    && apt-get install -y sudo \
    && rm -rf /var/lib/apt/lists/* \
    && echo "${USER_NAME:?} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers; \
    fi

# Return back to the project directory for the entrypoint
WORKDIR "${PROJECT_DIRECTORY}"
USER "${USER_NAME}"

# Prepare the directory for the VSCode server extensions and bash history
RUN mkdir -p "${USER_DIRECTORY}/.vscode-server/extensions" \
    "${USER_DIRECTORY}/.vscode-server-insiders/extensions" \
    && chown -R "${USER_NAME}:${USER_NAME}" \
    "${USER_DIRECTORY}/.vscode-server" \
    "${USER_DIRECTORY}/.vscode-server-insiders" && \
    SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=${USER_DIRECTORY}/.commandhistory/.bash_history" \
    && mkdir "${USER_DIRECTORY}/.commandhistory/" \
    && touch "${USER_DIRECTORY}/.commandhistory/.bash_history" \
    && chown -R "${USER_NAME}:${USER_NAME}" "${USER_DIRECTORY}/.commandhistory/" \
    && echo "$SNIPPET" >> "/home/$USER_NAME/.bashrc"

# Setup the pre-commit hooks so that they are run before each commit
RUN pre-commit install --install-hooks
