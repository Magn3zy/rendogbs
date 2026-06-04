FROM alpine:latest

RUN apk add --no-cache python3

COPY rendogbs_run.sh /home/rendogbs/
RUN chmod +x /home/rendogbs/rendogbs_run.sh
COPY rendogbs_pipeline.py /home/rendogbs/
COPY rendogbs_plots.py /home/rendogbs/
COPY endonucleases.py /home/rendogbs/

ENTRYPOINT ["/home/rendogbs/rendogbs_run.sh"]

