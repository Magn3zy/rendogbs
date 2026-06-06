.PHONY: all
all: .docker-built

.docker-built: Dockerfile src/endonucleases.py		\
	src/rendogbs_pipeline.py src/rendogbs_plots.py	\
	src/rendogbs_run.sh
	docker build -t rendogbs-v1 .
	touch $@

.PHONY: clean
clean:
	rm -f .docker-built
	docker image rm rendogbs-v1
