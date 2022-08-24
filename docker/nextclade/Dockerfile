# TAG dnastack/monkeypox-processing-pipeline/nextclade:2.4.0
FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@danstack.com>

RUN apt-get -qq update && \
  apt-get -qq install \
  wget \
  jq

ENV NEXTCLADE_VERSION 2.4.0
RUN mkdir -p /opt/nextclade && \
  wget -O /opt/nextclade/nextclade https://github.com/nextstrain/nextclade/releases/download/${NEXTCLADE_VERSION}/nextclade-x86_64-unknown-linux-gnu && \
 chmod +x /opt/nextclade/nextclade
RUN ln -s /opt/nextclade/nextclade /usr/local/bin

ENV MONKEYPOX_DATASET /data/monkeypox
RUN mkdir -p ${MONKEYPOX_DATASET} && \
  nextclade dataset get --name 'MPXV' --output-dir ${MONKEYPOX_DATASET}
