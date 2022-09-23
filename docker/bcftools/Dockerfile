# TAG dnastack/bcftools:1.16
FROM ubuntu:xenial

MAINTAINER Karen Fang <karen@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		build-essential \
		wget \
		tar \
		autoconf \
		bzip2 \
		lbzip2 \
		libbz2-dev \
  		libcurl4-openssl-dev \
  		liblzma-dev \
  		ncurses-dev \
 		libncurses5-dev \
		zlib1g-dev \
		libbz2-dev \
		jq \
		curl && \
	rm -rf /var/lib/apt/lists/*

# bcftools
ENV BCFTOOLS_VERSION 1.16
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
	tar -xvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt && \
	rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN	cd /opt/bcftools-${BCFTOOLS_VERSION} && \
	autoheader && \
	autoconf && \
	./configure && \
	make && \
	make install
ENV PATH $PATH:/opt/bcftools-${BCFTOOLS_VERSION}