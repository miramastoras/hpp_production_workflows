FROM tpesout/hpp_base:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

WORKDIR /opt

RUN git clone https://github.com/lh3/yak && \
	cd yak && \
	make
ENV PATH="/opt/yak:$PATH"

WORKDIR /data
