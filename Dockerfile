# Dockerfile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Docker image for rendogbs. Two-stage build with compiled binaries
# being created during build stage and included in the final image
# without the development environment.

################################################################
# Build Stage

FROM alpine:latest AS build

RUN apk add --no-cache rust cargo

COPY ./src/main.rs /home/build/
COPY ./src/Cargo.toml /home/build/

RUN cd /home/build ; cargo build --release

################################################################
# Runtime Stage

FROM alpine:latest

RUN apk add --no-cache python3 py3-numpy py3-matplotlib

COPY ./src/rendogbs_user.sh /home/rendogbs/
RUN chmod +x /home/rendogbs/rendogbs_user.sh

COPY ./src/rendogbs_run.sh /home/rendogbs/

COPY ./src/rendogbs_pipeline.py /home/rendogbs/
COPY ./src/rendogbs_plots.py /home/rendogbs/
COPY ./src/endonucleases.py /home/rendogbs/

COPY --from=build /home/build/target/release/rendogbs_finder /home/rendogbs/

ENTRYPOINT ["/home/rendogbs/rendogbs_user.sh"]
