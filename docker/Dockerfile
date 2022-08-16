FROM ubuntu:22.04

LABEL stage=builder

LABEL maintainer="tjakobi@arizona.edu"

ARG CIRCTOOLS_VERSION

ARG BUILD_DATE

LABEL org.label-schema.build-date=$BUILD_DATE

# update all repos, install minimal start packages

RUN sed -i 's/http:\/\/archive\.ubuntu\.com/http:\/\/mirror\.arizona\.edu/' /etc/apt/sources.list; sed -i 's/http:\/\/security\.ubuntu\.com/http:\/\/mirror\.arizona\.edu/' /etc/apt/sources.list; apt-get update; apt-get install --no-install-recommends wget git gpg ca-certificates -y

RUN useradd -ms /bin/bash circtools

# Set up repositories first

# RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | tee /etc/apt/trusted.gpg.d/R.gpg  >/dev/null; wget -qO- "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0xba6932366a755776" | gpg --dearmor | tee /etc/apt/trusted.gpg.d/python.gpg  >/dev/null; echo "deb [signed-by=/etc/apt/trusted.gpg.d/R.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list; echo "deb [signed-by=/etc/apt/trusted.gpg.d/python.gpg] https://ppa.launchpadcontent.net/deadsnakes/ppa/ubuntu jammy main" | tee -a /etc/apt/sources.list.d/python_ppa.list; mkdir /root/.R/;  wget https://raw.githubusercontent.com/jakobilab/circtools/devel/docker/Makevars -O /root/.R/Makevars; apt-get update -y

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | tee /etc/apt/trusted.gpg.d/R.gpg  >/dev/null; echo "deb [signed-by=/etc/apt/trusted.gpg.d/R.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list;  mkdir /root/.R/;  wget https://raw.githubusercontent.com/jakobilab/circtools/devel/docker/Makevars -O /root/.R/Makevars; apt-get update -y

# install R and Python

RUN apt-get update; DEBIAN_FRONTEND=noninteractive TZ=America/Phoenix apt-get install --no-install-recommends r-base python3 python3-dev make g++ gfortran libpng-dev zlib1g-dev libbz2-dev libjpeg-turbo8-dev libcurl4-openssl-dev libxml2-dev libblas-dev liblzma-dev libfontconfig1-dev liblapack-dev libssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libtiff5-dev libjpeg-dev -y;

# clone circtools
# install circtools deps (i.e. bedtools)
# install R packages
# install pip
# install circtools

RUN git clone -b devel https://github.com/jakobilab/circtools.git ; ls -la circtools/; bash /circtools/scripts/install_external.sh; Rscript /circtools/scripts/install_R_dependencies.R;

RUN wget https://bootstrap.pypa.io/get-pip.py; python3 get-pip.py; rm get-pip.py; python3 -m pip install -U setuptools; python3 -m pip install circtools/; pip cache purge; apt-get purge python3-dev; apt-get autoremove; apt-get autoclean; rm -rf /var/lib/apt/lists/*; rm circtools/ -rf;

USER circtools

# define entrypoint

ENTRYPOINT ["circtools"]



