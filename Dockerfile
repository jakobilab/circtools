FROM ubuntu:24.04

LABEL stage=builder
LABEL maintainer="tjakobi@arizona.edu"

ARG MAKEFLAGS="-j4"
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Phoenix
ENV PATH="/circtools/bin:$PATH"

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y wget gpg ca-certificates && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
    > /etc/apt/sources.list.d/r-project.list

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    git  \
    make \
    cmake \
    bzip2  \
    rsync  \
    g++  \
    gfortran \
    r-base  \
    pandoc \
    python3 \
    python3-dev  \
    python3-pip  \
    python3-venv \
    libpng-dev \
    zlib1g-dev  \
    libbz2-dev  \
    libjpeg-turbo8-dev \
    libopenblas-dev  \
    libcurl4-openssl-dev  \
    libxml2-dev  \
    libblas-dev \
    liblzma-dev \
    libgit2-dev  \
    libfontconfig1-dev  \
    liblapack-dev  \
    libssl-dev \
    libharfbuzz-dev  \
    uuid-dev  \
    libmariadb-dev-compat libfribidi-dev \
    libfreetype6-dev  \
    libtiff5-dev  \
    libncurses-dev \
    libjpeg-dev \
    libgmp-dev  \
    libmpfr-dev  \
    libnlopt-dev  \
    libuv1-dev \ 
    nano \
    curl && \
    useradd -ms /bin/bash circtools && \
    mkdir -p /root/.R && \
    mkdir /build/

COPY Makevars /root/.R/Makevars
ADD . /build/circtools/

RUN python3 -m venv /circtools && \
    . /circtools/bin/activate && \
    pip install --upgrade pip setuptools wheel && \
    pip install psutil && \
    pip install /build/circtools/ --verbose && \
    pip cache purge && \
    circtools_install_R_dependencies /build/circtools


RUN . /circtools/bin/activate && \
    circtools_install_R_dependencies /build/circtools && \
    Rscript -e 'library(primex); cat("primex OK\n")'

RUN cd /build && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && make -j4 && cp bin/* /usr/local/bin/ && \
    cd /build && \
    git clone --depth=1 https://github.com/icebert/pblat.git && \
    cd pblat && make && cp pblat /usr/local/bin/ && \
    cd /build && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar xvf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && make && make install && \
    apt-get purge python3-dev -y && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /build/ /var/lib/apt/lists/*

ADD docker_path_wrapper.py /usr/local/bin/
RUN chmod +x /usr/local/bin/docker_path_wrapper.py && mkdir /host_os/ && mkdir /host_rel/

LABEL org.opencontainers.image.description="Official circtools Docker image"
ENTRYPOINT ["docker_path_wrapper.py"]