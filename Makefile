#
# Makefile
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Simple build system for rendogbs.
#

# Configuration
IMGNAME=rendogbs-v2

################################################################
# User Targets

# By default, build the docker image
.PHONY: default
default: docker

# For manually starting the docker image build
.PHONY: docker
docker: containers/.docker-rendogbs-built

# Build everything at once
.PHONY: all
all: docker singularity

# Build only the singularity image (actually it builds the Docker
# image as well as it is built from it).
.PHONY: singularity
singularity: $(IMGNAME).sif

# Manually build only the bedtools Docker image
.PHONY: docker-bedtools
docker-bedtools: containers/.docker-bedtools-built

# Manually build only the rust binaries Docker image
.PHONY: docker-rust
docker-rust: containers/.docker-rust-built

################################################################
# Dependencies

# Interpreted scripts used in the pipeline(s)
SCRIPTS=src/endonucleases.py src/rendogbs_pipeline.py		\
	src/rendogbs_plots.py src/combination_processing.py	\
	src/enzymes.csv

# Wrappers needed for running in docker - the outer wrapper
# "rendogbs.sh" is NOT a dependency of the image!
WRAPPERS=src/rendogbs_run.sh src/rendogbs_user.sh

# Compiled programs (need to be compiled during build stage and then
# included in the final image only in binary form - see Dockerfile for
# details)
PROGRAMS=src/rust/cut_merge/cut_merge.rs				\
		src/rust/cut_merge/Cargo.toml				\
		src/rust/fragment_generation/fragment_generation.rs	\
		src/rust/fragment_generation/Cargo.toml			\
		src/rust/rendogbs_finder/rendogbs_finder.rs		\
		src/rust/rendogbs_finder/Cargo.toml

################################################################
# Clean Targets

# A bit crude, but works
.PHONY: clean
clean:
	rm -f containers/.docker-rendogbs-built $(IMGNAME).sif
	docker image rm $(IMGNAME) || true

# Like clean but also removes the compiled bedtools docker image
.PHONY: distclean
distclean: clean
	rm -f containers/.docker-bedtools-built
	rm -f containers/.docker-rust-built
	docker image rm rendogbs-bedtools || true
	docker image rm rendogbs-rust || true

################################################################
# Singularity Image

$(IMGNAME).sif: containers/.docker-rendogbs-built
	singularity build -F $@ docker-daemon://$(IMGNAME):latest

################################################################
# Main Docker Image

# A hidden file representing successful docker image build
containers/.docker-rendogbs-built: containers/$(IMGNAME)/Dockerfile	\
		$(SCRIPTS) containers/.docker-bedtools-built		\
		$(WRAPPERS) containers/.docker-rust-built
	docker build -f $< -t $(IMGNAME) .
	touch $@

################################################################
# Build-Support Docker Images

# A hidden file representing successful bedtools docker image build
containers/.docker-bedtools-built: containers/rendogbs-bedtools/Dockerfile
	docker build -f $< -t rendogbs-bedtools .
	touch $@

# A hidden file representing successful rust docker image build
containers/.docker-rust-built: containers/rendogbs-rust/Dockerfile	\
		$(PROGRAMS)
	docker build -f $< -t rendogbs-rust .
	touch $@
