FROM alpine:latest

RUN apk add --no-cache python3

COPY rendogbs2.sh /home/rendogbs/
RUN chmod +x /home/rendogbs/rendogbs2.sh

ENTRYPOINT ["/home/rendogbs/rendogbs2.sh"]

