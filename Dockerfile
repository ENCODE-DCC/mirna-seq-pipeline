# Dockerfile for ENCODE-DCC mirna-seq-pipeline
FROM ubuntu:16.04
MAINTAINER Otto Jolanki

RUN apt-get update && apt-get install -y \
    python3-dev \
    python3-pip \
    wget \
    libkrb5-3
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install STAR/Samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install && rm ../zlib-1.2.11.tar.gz

# Install STAR 2.5.1b
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz
RUN cd STAR-2.5.1b && make STAR && rm ../2.5.1b.tar.gz
ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

# Install wigToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && chmod +x wigToBigWig

# Install cutadapt 1.7.1
RUN pip3 install cutadapt==1.7.1

RUN mkdir -p mirna-seq-pipeline/src
COPY /src mirna-seq-pipeline/src
ENV PATH="/software/mirna-seq-pipeline/src:${PATH}"

ENTRYPOINT ["/bin/bash", "-c"]