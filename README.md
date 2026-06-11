rendogbs pipeline
=================

Requirements
------------

- POSIX.1 shell
- make
- docker

Building / Installation
-----------------------

To build the Docker image locally, just type:

```sh
make
```

The build process builds a Docker image named `rendogbs-v1` locally as
the default make target is `docker`. The image is built from supplied
`Dockerfile` in the root of this repository.

To build a singularity image, use the `singularity` target:

```sh
make singularity
```

This creates the `rendogbs-v1.sif` file.

It is possible to build both Docker and singularity image by using the
`all` target, however as the singularity image is built from the
Docker image, it basically performs the same operations as the
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

To run the Docker image a wrapper script is provided:

```sh
sh rendogbs.sh [ARGS]
```

Or alternatively:

```sh
sh rendogbs.sh --docker [ARGS]
```

The wrapper script uses Docker by default. It is possible to make it
run the singularity container `rendogbs-v1.sif` in current
directory. As the container is built in the root of this repository,
the following suffices:

```sh
sh rendogbs.sh --singularity [ARGS]
```

Running the wrapper script with the `--help` argument will print usage
instructions:

```sh
sh rendogbs.sh --help
```

Make sure you read the usage instructions and understand all the
options. Although the pipeline makes best effort to do the right
thing(TM) and report any issues as specifically as possible, it is
important to understand how it is designed and which resources it
leverages. See the minimal example below for a start as well.

Running Examples
----------------

For all the examples we assume a `run2` subdirectory in the current
directory will be used as a workdir and the `data` subdirectory will
contain a reference file `GCF_000001735.3_TAIR10_genomic.fna.gz`. All
the examples will run 5 parallel threads and use size 200-400.

### Running the Docker Image

Using the wrapper script, it is pretty straightforward:

```sh
sh rendogbs.sh --ref ./data/GCF_000001735.3_TAIR10_genomic.fna.gz --workdir ./run2 --parallel 5 --size 200-400
```

### Running the Singularity Image

The wrapper script handles all of the differences Singularity brings:

```sh
sh rendogbs.sh --singularity --ref ./data/GCF_000001735.3_TAIR10_genomic.fna.gz --workdir ./run2 --parallel 5 --size 200-400
```

### Running the Docker Image Manually

It is possible, but **strongly** discouraged, to run the Docker image
manually.

For Docker both the working directory and any and all files used by
the pipeline must be mapped as overlay volumes using the `-v`
option. It is advisable to bind the source files directly into the
runtime directory `/home/rendogbs` inside the container. When the
container terminates it should be removed by the `--rm` option. A full
example equivalent to the above examples would be:

```sh
docker run \
  -v ./GCF_000001735.3_TAIR10_genomic.fna.gz:/home/rendogbs/GCF_000001735.3_TAIR10_genomic.fna.gz \
  -v ./run2:/home/rendogbs/run2 \
  --rm \
  rendogbs-v1 \
  --ref ./GCF_000001735.3_TAIR10_genomic.fna.gz \
  --workdir ./run2 \
  --parallel 5 \
  --size 200-400
```

### Running the Singularity Container Manually

It is possible, but **strongly** discouraged, to run the Singularity
image manually.

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

Makefile
containers/rendogbs-v1/Dockerfile - Docker image recipy
containers/rendogbs-bedtools/Dockerfile - ...
rendogbs.sh - outer wrapper script
src/ - sources
src/endonucleases.py - default data source
src/rendogbs_pipeline.py - basic pipeline
src/rendogbs_plots.py - plots pipeline
src/rendogbs_run.sh - inner runtime wrapper script
src/rendogbs_user.sh - inner permission wrapper script

### Outer Wrapper Script: rendogbs.sh

TODO IMGNAME

This script validates the presence of mandatory arguments on the
command-line and creates appropriate Docker volume mappings for any
files and/or directories the pipeline needs.

Then it runs the Docker image `rendogbs-v1` and passes all the
collected arguments to its entrypoint which is the inner wrapper
script.

### Build System: Makefile

TODO IMGNAME

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

For Docker it should receive `LUID` and `LGID` (local user ID and
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

