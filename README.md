rendogbs pipeline
=================

Requirements
------------

- POSIX.1 shell
- make
- docker

Building / Installation
-----------------------

To build the docker image locally, just type:

```sh
make
```

The build process builds a docker image named `rendogbs-v1`
locally. The image is built from supplied `Dockerfile` in the root of
this repository.

Running
-------

```sh
sh rendogbs.sh [ARGS]
```

Architecture
------------

Dockerfile - docker image recipy
rendogbs.sh - outer wrapper script
src/ - sources
src/endonucleases.py - data
src/rendogbs_pipeline.py - basic pipeline
src/rendogbs_plots.py - plots pipeline
rendogbs_run.sh - inner wrapper script

### Outer Wrapper Script: rendogbs.sh

This script validates the presence of mandatory arguments on the
command-line and creates appropriate docker volume mappings for any
files and/or directories the pipeline needs.

Then it runs the docker image `rendogbs-v1` and passes all the
collected arguments to its entrypoint which is the inner wrapper
script.

### Docker Image: Dockerfile

Based on alpine Linux image it adds necessary Python version and
libraries used (numpy, matplotlib). It copies all the scripts to the
runtime directory inside the image (`/home/rendogbs`) and ensures the
inner wrapper script is used as image entrypoint.
