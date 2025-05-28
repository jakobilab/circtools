
FROM ubuntu:22.04

LABEL stage=builder
LABEL maintainer="tjakobi@arizona.edu"

ARG MAKEFLAGS="-j4"
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
    python3-venv \
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
    libgit2-dev \
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

# Download and install BEDTools
RUN cd /build/ && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make -j4 && \
    cp bin/* /usr/local/bin/ && \
    cd /build/ && \
    git clone --depth=1 https://github.com/icebert/pblat.git && \
    cd pblat && \
    make && \
    cp pblat /usr/local/bin/ &&\
    cd /build/ && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar xvf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    make && \
    make install &&\
    cd /build/ && \
    git clone --depth=1 https://github.com/ucscGenomeBrowser/kent.git && \
    cd kent/src/ && \
    make userApps && \
    cp ~/bin/`uname -m`/liftOver /usr/local/bin

# Download and install circtools
ADD . /build/circtools/

#RUN cd /build/ && \
#    python3 -m pip install -U setuptools numpy --break-system-packages && \
#    python3 -m pip install circtools/ --break-system-packages && \
#    cd /build/ && \
#    Rscript circtools/circtools/scripts/install_R_dependencies.R circtools/circtools/ &&\
#    pip install nanofilt --break-system-packages -v &&\
#    pip cache purge && \
#    apt-get purge python3-dev -y && \
#    apt-get autoremove -y && \
#    apt-get autoclean -y && \
#    rm /build/ /var/lib/apt/lists/ -rf


RUN python3 -m venv /circtools && \
    cd /build/ && \
    . /circtools/bin/activate && \
    pip install --upgrade pip setuptools wheel && \
    pip install Cython && \
    pip install numpy \
                primer3-py \
                "biopython>=1.71" && \
    pip install circtools/ --verbose && \
    cd /build/ && \
    Rscript circtools/circtools/scripts/install_R_dependencies.R circtools/circtools/ && \
    R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); \
          BiocManager::install(c('ballgown', 'ggbio', 'edgeR', 'GenomicRanges', 'GenomicFeatures'))" && \
    R -e "install.packages(c( \
      'ggplot2', 'ggrepel', 'plyr', 'ggfortify', 'openxlsx', \
      'formattable', 'kableExtra', 'dplyr', 'RColorBrewer', \
      'colortools', 'data.table', 'reshape2', 'gridExtra' \
    ), repos='https://cloud.r-project.org')" && \
    pip install nanofilt -v && \
    pip cache purge && \
    apt-get purge python3-dev -y && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /build/ /var/lib/apt/lists/*



# add script to bend absolute path names for circtools inside docker
ADD docker_path_wrapper.py /usr/local/bin/

RUN mkdir /host_os/

LABEL org.opencontainers.image.description="Official circtools Docker image"

# define entrypoint
ENTRYPOINT ["docker_path_wrapper.py"]
