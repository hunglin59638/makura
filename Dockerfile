FROM ubuntu:20.04

LABEL base.image="ubuntu"
LABEL software="makura"
LABEL software.version="1.1.0"
LABEL description="Makura: NCBI Genome downloader"
LABEL website="https://github.com/hunglin59638/makura"
LABEL maintainer="Hung-Lin Chen"
LABEL maintainer.email="hunglin59638@gmail.com"

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    rsync

RUN mkdir -p /opt/makura 
WORKDIR /opt/makura
ADD . .
RUN python3 setup.py install && \
    python3 makura-runner.py update --assembly-source refseq && \
    python3 makura-runner.py update --assembly-source genbank 
