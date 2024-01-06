FROM ubuntu:22.04

### PACKAGES ###

RUN apt-get update && apt-get install -y \
    libalberta-dev \
    libmetis-dev \
    libopenmpi-dev \
    libparmetis-dev \
    libsuitesparse-dev \
    git \
    g++ \
    cmake \
    pkg-config

### DUNE / AMDIS ###

RUN mkdir /dune && \
    cd /dune && \
    git clone https://gitlab.dune-project.org/core/dune-common && \
    git clone https://gitlab.dune-project.org/core/dune-geometry && \
    git clone https://gitlab.dune-project.org/core/dune-grid && \
    git clone https://gitlab.dune-project.org/core/dune-istl && \
    git clone https://gitlab.dune-project.org/core/dune-localfunctions && \
    git clone https://gitlab.dune-project.org/staging/dune-functions && \
    git clone https://gitlab.dune-project.org/staging/dune-typetree && \
    git clone https://gitlab.dune-project.org/staging/dune-uggrid && \
    git clone https://gitlab.com/amdis/amdis.git amdis

RUN cd /dune && \
    export CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo" && \
    export MAKE_FLAGS="-j4" && \
    dune-common/bin/dunecontrol all

### PROJECT DIR ###

RUN mkdir /bio_condensates && cd /bio_condensates
WORKDIR /bio_condensates

COPY . .

ENTRYPOINT ["scripts/run.sh"]
