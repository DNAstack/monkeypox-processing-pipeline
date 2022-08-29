# TAG dnastack/bwa_samtools:0.0.1
FROM ubuntu:xenial

MAINTAINER Karen Fang <karen@dnastack.com>

RUN apt-get -qq update && \
	apt-get -qq install \
		build-essential \
		wget \
		tar \
		gzip \
		tabix \
		lbzip2 \
		make \
		autoconf \
		libncurses5-dev \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev && \
	rm -rf /var/lib/apt/lists/*

# bwa
ENV BWA_VERSION 0.7.17
RUN wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 && \
	tar -xvf bwa-${BWA_VERSION}.tar.bz2 --directory /opt && \
	rm bwa-${BWA_VERSION}.tar.bz2
RUN	cd /opt/bwa-${BWA_VERSION} && \
	make
ENV PATH $PATH:/opt/bwa-${BWA_VERSION}

# samtools
ENV SAMTOOLS_VERSION 1.15
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar -xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt && \
	rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN	cd /opt/samtools-${SAMTOOLS_VERSION} && \
	autoheader && \
	autoconf -Wno-syntax && \
	./configure && \
	make && \
	make install
