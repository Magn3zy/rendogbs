FROM alpine:latest

RUN apk add --no-cache python3

COPY rendogbs2.sh /home/rendogbs/
RUN chmod +x /home/rendogbs/rendogbs2.sh
COPY rendogbs_pipeline.py /home/rendogbs/
COPY rendogbs_plots.py /home/rendogbs/
COPY endonucleases.py /home/rendogbs/

ENTRYPOINT ["/home/rendogbs/rendogbs2.sh"]

