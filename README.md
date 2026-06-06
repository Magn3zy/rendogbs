rendogbs pipeline
=================

Requirements
------------

- docker

Building / Installation
-----------------------

```sh
docker build -t rendogbs01 .
```

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

