FROM ubuntu:18.04

LABEL \
  version="1.12" \
  description="bcftools image for use in Workflows"

RUN apt-get update && apt-get install -y \
  bzip2 \
  g++ \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev

ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV BCFTOOLS_VERSION=1.12

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 && \
  tar --bzip2 -xf bcftools-$BCFTOOLS_VERSION.tar.bz2

WORKDIR /tmp/bcftools-$BCFTOOLS_VERSION
RUN make prefix=$BCFTOOLS_INSTALL_DIR && \
  make prefix=$BCFTOOLS_INSTALL_DIR install

WORKDIR /
RUN ln -s $BCFTOOLS_INSTALL_DIR/bin/bcftools /usr/bin/bcftools && \
  rm -rf /tmp/bcftools-$BCFTOOLS_VERSION

# put in path
ENV PATH=${PATH}:/usr/bin/bcftools

RUN apt-get update && \
	apt-get install -y vcftools


ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev awscli jq pigz \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 --no-cache-dir install --upgrade pip \
  && rm -rf /var/lib/apt/lists/* \
  && python -m pip install --upgrade awscli

RUN python -m pip install boto3 pandas pysam plotnine

WORKDIR /opt/

RUN apt-get update && \
    apt-get -y install time git wget && \
    apt-get clean && \
    apt-get purge

## Get script

RUN git clone https://github.com/miramastoras/useful_scripts.git

ENV PATH=/opt/useful_scripts/rmdups_vcf.py:$PATH
