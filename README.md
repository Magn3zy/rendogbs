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

The build process builds a Docker image named `rendogbs-v2` locally as
the default make target is `docker`. The image is built from supplied
`Dockerfile` in the root of this repository.

To build a singularity image, use the `singularity` target:

```sh
make singularity
```

This creates the `rendogbs-v2.sif` file.

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
run the singularity container `rendogbs-v2.sif` in current
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
  rendogbs-v2 \
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
  ./rendogbs-v2.sif \
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
  ./rendogbs-v2.sif \
  --ref ./GCF_000001735.3_TAIR10_genomic.fna.gz \
  --workdir ./run2 \
  --parallel 5 \
  --size 200-400
```

These two options are equivalent.


Argument guide
--------------
```
REQUIRED
--ref            <file>   Reference FASTA (.fa / .fasta / .fa.gz)
--workdir        <dir>    Working directory
--parallel       <int>    Combinations processed in parallel per batch, don't use more threads than combinations
--size           <range>  Fragment size window used in library, e.g. 200-400 (both inclusive).
                          Standard 100 bp bins (0-99 .. 900-999 + >=1000) are
                          always reported; filtered.csv retains only fragments
                          within --size-range.
--chroms         <int>    Longest N contigs treated as chromosomes in
                          per-chromosome plots (default: 10)

COMBINATIONS (fast/custom)
--combinations            fast      (deafult most used combinations from literature)
                                    you don't need to specify this argument
--combinations            custom    (specify combinations yourself)
--combinations-file       <file>    (specify combinations you csv file) 
                                    csv formating: header line - enzyme_a,enzyme_b, refer to enzyme list

ANNOTATION (optional)
--annotation     <file>   GFF3/GFF/GTF gene annotation
--te             <file>   RepeatMasker .out TE annotation

PLOT OPTIONS (optional)
--dpi            <int>    Figure resolution in DPI (default: 300)
--skip-plots              Run pipeline only, skip plot generation

HELP
-h, --help                Show help (shows commented file structure used)
--enzymes                 Show available enzymes in this pipeline (171)
--fast-combinations       Show fast combinations used in this pipeline
```
	  

Architecture
------------

The whole pipeline is packaged as a container which includes all the
scripts, programs, and their dependencies. There are two supported
container platforms: Docker and Apptainer (formerly Singularity).

As the pipeline needs access to local files, proper file and directory
mappings need to be provided for either type of container. The Outer
Wrapper Script runs the container with mappings after it performs
preliminary argument validation.

The entry point of the container is the Inner User Wrapper Script
which is needed under Docker containerization to setup the the user
and group under which the actual pipeline runs to match the outer
environment. With Apptainer this layer does not perform any
adjustments as the environment is already correctly set up. Then the
Inner Pipeline Wrapper Script is run.

The Inner Pipeline Wrapper Script validates all the arguments and
assigns them to appropriate argument sets of individual steps. Then it
performs all the steps and measures their running times. If any of the
steps fails, the whole pipeline fails immediately and such information
is reported to the user.

Upon successful pipeline run timing information is shown and the
container exits successfully.

### Build System: Makefile, Dockerfiles and SIF Conversion

Standard Makefile with target dependencies is used for implementing
the build process. The default target builds the Docker image and a
separate `sif` target builds the Apptainer image.

Internally there are multiple interdependent targets that orchestrate
the whole process. For building the primary Docker image, two helper
images are built.

Firstly as the internal `docker-bedtools` Makefile target the
`rendogbs-bedtools` container is built using the
`containers/rendogbs-bedtools/Dockerfile` recipy. It builts the
`bedtools` binary using an older Alpine Linux 3.18 based image as the
version this project uses can only be built with gcc 12.

Secondly as the internal `docker-rust` Makefile target the
`rendogbs-rust` container is built which compiles all the Rust
binaries implemented for this project. It uses the latest Alpine Linux
image available as the latest Rust compiler is always needed.

Thirdly as the internal `docker` Makefile target the final Docker
image `rendogbs-v2` container is built from latest Alpine Linux
image. The bedtools binary and all the Rust programs are copied from
the former two containers without development files, include headers
and static libraries, reducing significantly the final image size. All
the scripts, binaries and support data are stored in `/home/rendogbs`
directory inside the image. The entry point is the
`/home/rendogbs/rendogbs_user.sh` shell script.

The outer wrapper script expects this `rendogbs-v2` image locally
available.

If a `sif` target is build it ensures the `rendogbs-v2` docker image
is locally available by dependeing on the `docker` internal target. It
checks for the `apptainer` (preferred) or `singularity` binary
availability and uses it to build the resulting SIF image using the
tool found.

### Outer Wrapper Script: rendogbs.sh

TODO IMGNAME

This script validates the presence of mandatory arguments on the
command-line and creates appropriate Docker volume mappings for any
files and/or directories the pipeline needs.

Then it runs the Docker image `rendogbs-v2` and passes all the
collected arguments to its entrypoint which is the inner wrapper
script.


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

