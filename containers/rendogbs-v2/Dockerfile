# Dockerfile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Docker image for rendogbs. Multi-stage build with compiled binaries
# being created during build stage and included in the final image
# without the development environment.
#
# Requirements:
# - rendogbs-bedtools
# - rendogbs-rust

################################################################
# Runtime Stage

FROM alpine:latest

RUN apk add --no-cache python3 py3-numpy py3-matplotlib py3-pandas

COPY src/rendogbs_user.sh /home/rendogbs/
RUN chmod +x /home/rendogbs/rendogbs_user.sh

COPY src/rendogbs_run.sh /home/rendogbs/

COPY src/combination_processing.py /home/rendogbs/
COPY src/fragment_generation.py /home/rendogbs/

COPY src/rendogbs_pipeline.py /home/rendogbs/
COPY src/rendogbs_plots.py /home/rendogbs/
COPY src/endonucleases.py /home/rendogbs/
COPY src/enzymes.csv /home/rendogbs/
COPY src/allowed_pairs.csv /home/rendogbs/

COPY --from=rendogbs-rust:latest /home/build/rendogbs_finder/target/release/rendogbs_finder /home/rendogbs/
COPY --from=rendogbs-rust:latest /home/build/cut_merge/target/release/cut_merge /home/rendogbs/
COPY --from=rendogbs-rust:latest /home/build/fragment_generation/target/release/fragment_generation /home/rendogbs/
COPY --from=rendogbs-bedtools:latest /home/build/bedtools2/bin/bedtools /usr/bin/bedtools

ENTRYPOINT ["/home/rendogbs/rendogbs_user.sh"]
