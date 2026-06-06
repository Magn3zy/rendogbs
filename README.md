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

The build process builds a docker image named `rendogbs-v1` locally as
the default make target is `docker`. The image is built from supplied
`Dockerfile` in the root of this repository.

To build a singularity image, use the `singularity` target:

```sh
make singularity
```

This creates the `rendogbs-v1.sif` file.

It is possible to build both Docker and singularity image by using the
`all` target, however as the singularity image is built from the
docker image, it basically performs the same operations as the
`singularity` target:

```sh
make all
```

To remove all built images and build artifacts, use:

```sh
make clean
```

Running
-------

To run the docker image a wrapper script is provided:

```sh
sh rendogbs.sh [ARGS]
```

Or alternatively:

```sh
sh rendogbs.sh --docker [ARGS]
```

The wrapper script uses docker by default. It is possible to make it
run the singularity container `rendogbs-v1.sif` in current
directory. As the container is built in the root of this repository,
the following suffices:

```sh
sh rendogbs.sh --singularity [ARGS]
```

Running the wrapper script with the `--help` argument will print usage
instructions.

Running Examples
----------------

For all the examples we assume a `run2` subdirectory in the current
directory will be used as a workdir and the `data` subdirectory will
contain a reference file `GCF_000001735.3_TAIR10_genomic.fna.gz`. All
the examples will run 5 parallel threads and use size 200-400.

### Running the Docker Image

```sh
```

### Running the Singularity Image

```sh
```

### Running the Docker Image Manually

It is possible, but *strongly* discouraged, to run the docker image
manually.

### Running the Singularity Container Manually

There are two options. The first option is to use the
`SINGULARITY_BIND` environment variable to bind files and directories
inside the container and run the container file directly as it is a
valid executable:

```sh
SINGULARITY_BIND=./run2:/home/rendogbs/run2,./data/GCF_000001735.3_TAIR10_genomic.fna.gz:/home/rendogbs/GCF_000001735.3_TAIR10_genomic.fna.gz \
  ./rendogbs-v1.sif \
    --ref ./GCF_000001735.3_TAIR10_genomic.fna.gz \
	--workdir ./run2 \
	--parallel 5 \
	--size 200-400 
```

As we can see there is no need to bind the whole `data` subdirectory
and it is easier to bind the reference file directly in the
`/home/rendogbs` runtime directory.

The other option is to use the `singularity run` command to provide
the bindings through the `-B` option:

```sh
singularity run \
  -B ./run2:/home/rendogbs/run2 \
  -B ./data/GCF_000001735.3_TAIR10_genomic.fna.gz:/home/rendogbs/GCF_000001735.3_TAIR10_genomic.fna.gz \
  ./rendogbs-v1.sif \
  --ref ./GCF_000001735.3_TAIR10_genomic.fna.gz \
  --workdir ./run2 \
  --parallel 5 \
  --size 200-400
```

These two options are equivalent.


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

### Inner User Wrapper Script: rendogbs_user.sh

This is a thin layer needed for typical Docker usage which ensures
that all the files are created with the current user as their
owner. As it is the entry point it handles setting up the environment
for both Docker and Singularity variangs.

For docker it should receive `LUID` and `LGID` (local user ID and
local group ID) environment variables and it will setup the inner
container environment to reflect these. The whole pipeline is then run
as user with the same UID and GID as those provided and all the files
have their ownership updated accordingly. When running manually under
Docker, the following arguments to `docker run` should always be
present:

```sh
docker run -e LUID=$(id -u) -e LGID=$(id -g) ...
```

For Singularity it receives the `SINGULARITY_CONTAINER` environment
variable and recognizes it is already running as the correct
unprivileged user.

### Inner Pipeline Wrapper Script: rendogbs_run.sh

TODO ...

