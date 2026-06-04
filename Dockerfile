FROM alpine:latest

RUN apk add --no-cache python3 py3-numpy py3-matplotlib

COPY ./src/rendogbs_run.sh /home/rendogbs/
RUN chmod +x /home/rendogbs/rendogbs_run.sh
COPY ./src/rendogbs_pipeline.py /home/rendogbs/
COPY ./src/rendogbs_plots.py /home/rendogbs/
COPY ./src/endonucleases.py /home/rendogbs/

ENTRYPOINT ["/home/rendogbs/rendogbs_run.sh"]

