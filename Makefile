#
# Makefile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Simple build system for rendogbs.
#

# By default, build the docker image
.PHONY: default
default: docker

# Build everything at once
.PHONY: all
all: docker singularity

# Interpreted scripts used in the pipeline(s)
SCRIPTS=src/endonucleases.py src/rendogbs_pipeline.py		\
	src/rendogbs_plots.py src/combination_processing.py	\
	enzymes.csv

# Wrappers needed for running in docker - the outer wrapper
# "rendogbs.sh" is NOT a dependency of the image!
WRAPPERS=src/rendogbs_run.sh src/rendogbs_user.sh

# Compiled programs (need to be compiled during build stage and then
# included in the final image only in binary form - see Dockerfile for
# details)
PROGRAMS=src/main.rs src/Cargo.toml

# For manually starting the docker image build
.PHONY: docker
docker: .docker-built

# A hidden file representing successful docker image build
.docker-built: containers/rendogbs-v1/Dockerfile $(SCRIPTS)	\
		.docker-bedtools-built $(WRAPPERS) $(PROGRAMS)
	docker build -f $< -t rendogbs-v1 .
	touch $@

# A bit crude, but works
.PHONY: clean
clean:
	rm -f .docker-built rendogbs-v1.sif
	docker image rm rendogbs-v1

# Like clean but also removes the compiled bedtools docker image
.PHONY: distclean
distclean: clean
	rm -f .docker-bedtools-built
	docker image rm rendogbs-bedtools

.PHONY: singularity
singularity: rendogbs-v1.sif

rendogbs-v1.sif: .docker-built
	singularity build -F rendogbs-v1.sif docker-daemon://rendogbs-v1:latest

# For manually starting the bedtools docker image build
.PHONY: docker-bedtools
docker-bedtools: .docker-bedtools-built

# A hidden file representing successful bedtools docker image build
.docker-bedtools-built: containers/rendogbs-bedtools/Dockerfile
	docker build -f $< -t rendogbs-bedtools .
	touch $@
