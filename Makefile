#
# Makefile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Simple build system for rendogbs.
#

# By default, build the docker image
.PHONY: all
all: docker

# Interpreted scripts used in the pipeline(s)
SCRIPTS=src/endonucleases.py src/rendogbs_pipeline.py	\
	src/rendogbs_plots.py

# Wrappers needed for running in docker - the outer wrapper
# "rendogbs.sh" is NOT a dependency of the image!
WRAPPERS=src/rendogbs_run.sh src/rendogbs_user.sh

# Compiled programs (need to be compiled during build stage and then
# included in the final image only in binary form - see Dockerfile for
# details)
PROGRAMS=src/program1.rs

.PHONY: docker
docker: .docker-built

# A hidden file representing successful docker image build
.docker-built: Dockerfile $(SCRIPTS) $(WRAPPERS) $(PROGRAMS)
	docker build -t rendogbs-v1 .
	touch $@

# A bit crude, but works
.PHONY: clean
clean:
	rm -f .docker-built rendogbs-v1.sif
	docker image rm rendogbs-v1

.PHONY: singularity
singularity: rendogbs-v1.sif

rendogbs-v1.sif: .docker-built
	singularity build rendogbs-v1.sif docker-daemon://rendogbs-v1:latest
