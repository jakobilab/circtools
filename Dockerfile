
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

ADD . /build/circtools/

# Install Git
RUN apt-get update && apt-get install -y git

# Step 1: Install Git
RUN apt-get update && apt-get install -y git

# Step 2: Install BiocManager
RUN Rscript -e "\
  cat('>>> Installing BiocManager...\n'); \
  if (!requireNamespace('BiocManager', quietly = TRUE)) \
    install.packages('BiocManager', repos='https://cloud.r-project.org'); \
  stopifnot(requireNamespace('BiocManager')); \
  cat('>>> BiocManager installed ✅\n')"


# Install math libraries needed by lme4, car, etc.
RUN apt-get update && apt-get install --no-install-recommends -y \
    libgmp-dev \
    libmpfr-dev \
    libnlopt-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*



# Step 3: Install CRAN deps for GGally, Hmisc, and ggstats
# Step A: Core packages
RUN Rscript -e "\
  install.packages(c('forcats', 'reshape', 'remotes', 'tibble', 'rlang', 'stringr'), \
                   repos='https://cloud.r-project.org')"

# Step B: Heavy deps (after system libs)
RUN Rscript -e "\
  install.packages(c('ggplot2', 'scales', 'survey', 'car', 'vcd'), \
                   repos='https://cloud.r-project.org')"



RUN Rscript -e "\
  cat('>>> Installing CRAN packages: Formula, latticeExtra, htmlTable, viridis...\n'); \
  install.packages(c('Formula', 'latticeExtra', 'htmlTable', 'viridis'), repos='https://cloud.r-project.org'); \
  stopifnot(all(c('Formula', 'latticeExtra', 'htmlTable', 'viridis') %in% rownames(installed.packages()))); \
  cat('>>> Core CRAN packages installed ✅\n')"

# Step 4: Install archived Hmisc (ensure data.table is present first)
RUN Rscript -e "\
  cat('>>> Installing data.table (needed by Hmisc)...\n'); \
  install.packages('data.table', repos='https://cloud.r-project.org'); \
  stopifnot('data.table' %in% rownames(installed.packages())); \
  cat('>>> data.table installed ✅\n')"

RUN Rscript -e "\
  cat('>>> Installing archived Hmisc...\n'); \
  install.packages('https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.6-0.tar.gz', \
                   repos=NULL, type='source'); \
  stopifnot('Hmisc' %in% rownames(installed.packages())); \
  cat('>>> Hmisc installed ✅\n')"


# Step: Install ggstats dependencies
RUN Rscript -e "\
  cat('>>> Installing ggstats dependencies: broom.helpers, dplyr, tidyr...\n'); \
  install.packages(c('broom.helpers', 'dplyr', 'tidyr'), repos = 'https://cloud.r-project.org'); \
  cat('>>> Dependencies installed ✅\n')"

# Step: Install archived ggstats 0.3.0
RUN Rscript -e "\
  cat('>>> Installing archived ggstats 0.3.0...\n'); \
  install.packages('https://cran.r-project.org/src/contrib/Archive/ggstats/ggstats_0.3.0.tar.gz', \
                   repos = NULL, type = 'source'); \
  stopifnot('ggstats' %in% rownames(installed.packages())); \
  cat('>>> ggstats 0.3.0 installed ✅\n')"


# Step 6: Install archived GGally
RUN Rscript -e "\
  cat('>>> Installing archived GGally...\n'); \
  install.packages('https://cran.r-project.org/src/contrib/Archive/GGally/GGally_2.1.2.tar.gz', repos=NULL, type='source'); \
  stopifnot('GGally' %in% rownames(installed.packages())); \
  cat('>>> GGally installed ✅\n')"

# Step 7: Install Bioconductor base packages for ggbio
RUN Rscript -e "\
  cat('>>> Installing Bioconductor base packages...\n'); \
  BiocManager::install(c('S4Vectors', 'IRanges', 'GenomicRanges', 'Biostrings', \
                         'BSgenome', 'GenomicFeatures', 'edgeR', 'ballgown'), \
                         ask=FALSE, update=FALSE); \
  stopifnot('IRanges' %in% rownames(installed.packages())); \
  cat('>>> Bioconductor base installed ✅\n')"

# Step 8: Install biovizBase
RUN Rscript -e "\
  cat('>>> Installing biovizBase...\n'); \
  BiocManager::install('biovizBase', ask=FALSE, update=FALSE); \
  stopifnot('biovizBase' %in% rownames(installed.packages())); \
  cat('>>> biovizBase installed ✅\n')"

# Step 9: Install ggbio and verify
RUN Rscript -e "\
  cat('>>> Installing ggbio...\n'); \
  BiocManager::install('ggbio', ask=FALSE, update=FALSE); \
  stopifnot('ggbio' %in% rownames(installed.packages())); \
  cat('>>> ggbio installed and verified ✅\n')"



# Step 6a: Install devtools (needed for GitHub installs)
RUN Rscript -e "\
  cat('>>> Installing devtools (for GitHub packages)...\n'); \
  install.packages('devtools', repos='https://cloud.r-project.org'); \
  cat('>>> devtools installed ✅\n')"

# Step 6: Install circtools extras from GitHub
RUN Rscript -e "\
  cat('>>> Installing circTest and primex from GitHub...\n'); \
  devtools::install_github('dieterich-lab/CircTest'); \
  devtools::install_github('dieterich-lab/primex'); \
  cat('>>> circTest and primex installed ✅\n')"



RUN apt-get update && apt-get install -y \
    libz-dev \
    libpthread-stubs0-dev \
    build-essential

# pblat
RUN cd /build && \
    git clone --depth=1 https://github.com/icebert/pblat.git && \
    cd pblat && \
    unset MAKEFLAGS && \
    make -j1 && \
    cp pblat /usr/local/bin/


# BEDTools
RUN cd /build && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/* /usr/local/bin/



# samtools
RUN cd /build && \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar xvf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    make && \
    make install

# UCSC liftOver binary (precompiled, avoids `make` failure on ARM/M1)
RUN apt-get update && apt-get install -y curl && \
    curl -LO https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod +x liftOver && \
    mv liftOver /usr/local/bin/


#primer docker build
RUN cd /build/ && \
    git clone https://github.com/primer3-org/primer3.git && \
    cd primer3/src && \
    make && \
    mkdir -p /usr/local/lib/R/site-library/primex/primer3 && \
    cp primer3_core /usr/local/lib/R/site-library/primex/primer3/primer3_core_Linux_64 && \
    chmod +x /usr/local/lib/R/site-library/primex/primer3/primer3_core_Linux_64


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


# Create circtools environment and install dependencies
RUN python3 -m venv /circtools && \
    . /circtools/bin/activate && \
    pip install --upgrade pip setuptools wheel && \
    pip install Cython numpy matplotlib "biopython>=1.71" primer3-py


# Install circtools
RUN . /circtools/bin/activate && \
    pip install /build/circtools/ --verbose

# Install StringTie from binary
RUN cd /build && \
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar -xvzf stringtie-2.2.1.Linux_x86_64.tar.gz && \
    cp stringtie-2.2.1.Linux_x86_64/stringtie /usr/local/bin/ && \
    chmod +x /usr/local/bin/stringtie









RUN R -e "install.packages(c( \
  'ggplot2', 'ggrepel', 'plyr', 'ggfortify', 'openxlsx', \
  'formattable', 'kableExtra', 'dplyr', 'RColorBrewer', \
  'colortools', 'data.table', 'reshape2', 'gridExtra' \
), repos='https://cloud.r-project.org')"

# Install remaining Python tools and cleanup
RUN . /circtools/bin/activate && \
    pip install nanofilt -v && \
    pip cache purge

RUN apt-get purge python3-dev -y && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /build/ /var/lib/apt/lists/*






# add script to bend absolute path names for circtools inside docker
ADD docker_path_wrapper.py /usr/local/bin/

RUN mkdir /host_os/

LABEL org.opencontainers.image.description="Official circtools Docker image"

# define entrypoint
ENTRYPOINT ["docker_path_wrapper.py"]
