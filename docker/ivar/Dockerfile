# TAG dnastack/monkeypox-processing-pipeline/ivar:1.3.1
FROM dnastack/monkeypox-processing-pipeline/bwa_samtools:0.0.1

MAINTAINER Heather Ward <heather@dnastack.com>

ENV HTSLIB_VERSION ${SAMTOOLS_VERSION}
RUN cd /opt/samtools-${SAMTOOLS_VERSION}/htslib-${HTSLIB_VERSION} && \
  ./configure && \
  make && \
  make install

ENV IVAR_VERSION 1.3.1
RUN wget https://github.com/andersen-lab/ivar/archive/refs/tags/v${IVAR_VERSION}.tar.gz && \
  tar -zxvf v${IVAR_VERSION}.tar.gz --directory /opt && \
  rm v${IVAR_VERSION}.tar.gz

RUN cd /opt/ivar-${IVAR_VERSION} && \
  ./autogen.sh && \
  ./configure && \
  make && \
  make install

ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/usr/local/lib"
