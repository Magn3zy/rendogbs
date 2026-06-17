# Dockerfile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Docker image for rendogbs-rust binaries. To be used for building the
# binaries before embedding in the final runtime image.
#

FROM alpine:latest AS build
RUN apk add --no-cache rust cargo

COPY ./src/rust/ /home/build/
RUN find /home/build

RUN cd /home/build/cut_merge && cargo build --release
RUN cd /home/build/fragment_generation && cargo build --release
RUN cd /home/build/rendogbs_finder && cargo build --release
