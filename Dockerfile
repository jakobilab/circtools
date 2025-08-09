FROM ubuntu:22.04

LABEL stage=builder
LABEL maintainer="tjakobi@arizona.edu"

ARG MAKEFLAGS="-j4"
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Phoenix
ENV PATH="/circtools/bin:$PATH"


RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    wget git gpg ca-certificates make bzip2 rsync g++ gfortran \
    r-base python3 python3-dev python3-pip python3-venv \
    libpng-dev zlib1g-dev libbz2-dev libjpeg-turbo8-dev \
    libopenblas-dev libcurl4-openssl-dev libxml2-dev libblas-dev \
    liblzma-dev libgit2-dev libfontconfig1-dev liblapack-dev libssl-dev \
    libharfbuzz-dev uuid-dev libmariadb-dev-compat libfribidi-dev \
    libfreetype6-dev libtiff5-dev libncurses-dev libjpeg-dev \
    libgmp-dev libmpfr-dev libnlopt-dev curl && \
    useradd -ms /bin/bash circtools && \
    mkdir -p /root/.R && \
    mkdir /build/

COPY Makevars /root/.R/Makevars
ADD . /build/circtools/


RUN python3 -m venv /circtools && \
    . /circtools/bin/activate && \
    pip install --upgrade pip setuptools wheel && \
    pip install /build/circtools/ --verbose && \
    pip cache purge && \
    Rscript /build/circtools/circtools/scripts/install_R_dependencies.R /build/circtools/

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
RUN chmod +x /usr/local/bin/docker_path_wrapper.py && mkdir /host_os/

LABEL org.opencontainers.image.description="Official circtools Docker image"
ENTRYPOINT ["docker_path_wrapper.py"]