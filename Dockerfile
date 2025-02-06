FROM ubuntu:22.04

LABEL stage=builder
LABEL maintainer="tjakobi@arizona.edu"

ARG MAKEFLAGS="-j1"
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Phoenix

# update all repos, install minimal start packages
RUN apt-get update && \
    apt-get install --no-install-recommends wget git gpg ca-certificates -y

RUN useradd -ms /bin/bash circtools

COPY Makevars /root/.R/Makevars

# Install all base packages for R and building
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive TZ=America/Phoenix apt-get install --no-install-recommends -y \
    r-base \
    python3  \
    python3-dev \
    python3-pip \
    make  \
    bzip2 \
    rsync \
    g++ \
    gfortran  \
    libpng-dev \
    zlib1g-dev \
    libbz2-dev \
    libjpeg-turbo8-dev \
    libopenblas-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libblas-dev \
    liblzma-dev \
    libfontconfig1-dev \
    liblapack-dev \
    libssl-dev \
    libharfbuzz-dev \
    uuid-dev \
    libmariadb-dev-compat \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libncurses-dev \
    libjpeg-dev

# make build dir
RUN mkdir /build/

## Download and install BEDTools
#RUN cd /build/ && \
#    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
#    tar zxvf bedtools-2.31.1.tar.gz && \
#    cd bedtools2 && \
#    make -j4 && \
#    cp bin/* /usr/local/bin/
#
## Download and install pblat
#RUN cd /build/ && \
#    git clone --depth=1 https://github.com/icebert/pblat.git && \
#    cd pblat && \
#    make && \
#    cp pblat /usr/local/bin/
#
## Download and install samtools
#RUN cd /build/ && \
#    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
#    tar xvf samtools-1.21.tar.bz2 && \
#    cd samtools-1.21 && \
#    make && \
#    make install

## Download and install circtools
#RUN cd /build/ && \
#    git clone --depth=1 https://github.com/jakobilab/circtools.git  && \
#    python3 -m pip install -U setuptools numpy --break-system-packages && \
#    python3 -m pip install circtools/ --break-system-packages

# Download and install circtools
ADD . /build/circtools/

#RUN cd /build/ && \
#    python3 -m pip install -U setuptools numpy --break-system-packages && \
#    python3 -m pip install circtools/ --break-system-packages
#
## Download and install circtools R deps
#RUN cd /build/ && \
#    Rscript circtools/circtools/scripts/install_R_dependencies.R circtools/circtools/
#
#RUN pip install nanofilt --break-system-packages -v

RUN uname -m

RUN cd /build/ && \
    git clone --depth=1 https://github.com/ucscGenomeBrowser/kent.git && \
    cd kent/src/ && \
    make userApps && \
    cp ~/bin/`uname -m`/liftOver /usr/local/bin

# Clean up to save space
RUN pip cache purge && \
    apt-get purge python3-dev -y && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm /build/ /var/lib/apt/lists/ -rf

# add script to bend absolute path names for circtools inside docker
ADD docker_path_wrapper.py /usr/local/bin/


RUN mkdir /host_os/

LABEL org.opencontainers.image.description="Official circtools Docker image"

# define entrypoint
ENTRYPOINT ["docker_path_wrapper.py"]
