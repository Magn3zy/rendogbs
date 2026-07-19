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

To build an Apptainer (formerly Singularity) image, use the `sif`
target:

```sh
make sif
```

This creates the `rendogbs-v2.sif` file.

It is possible to build both Docker and singularity image by using the
`all` target, however as the singularity image is built from the
Docker image, it basically performs the same operations as the
`sif` target:

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
run the Apptainer container `rendogbs-v2.sif` in current directory. As
the container is built in the root of this repository, the following
suffices:

```sh
sh rendogbs.sh --apptainer [ARGS]
```

For backwards compatibility, the following works as well:

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

### Running the Apptainer (formerly Singularity) Image

The wrapper script handles all of the differences Apptainer brings:

```sh
sh rendogbs.sh --apptainer --ref ./data/GCF_000001735.3_TAIR10_genomic.fna.gz --workdir ./run2 --parallel 5 --size 200-400
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

### Running the Apptainer Container Manually

It is possible, but **strongly** discouraged, to run the Apptainer
image manually. The following guide uses the legacy `singularity`
binary which works both in Apptainer and original Singularity
distribution.

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
	  