# Dockerfile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Docker image for compiling bedtools binary for rendogbs. This is a
# reusable build stage.


################################################################
# Build Stage - bedtools

FROM alpine:3.18
RUN apk add bash build-base zlib-dev bzip2-dev xz-dev python3

RUN mkdir /home/build

RUN cd /home/build && wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && tar zxf bedtools-2.31.1.tar.gz && cd bedtools2 && make
