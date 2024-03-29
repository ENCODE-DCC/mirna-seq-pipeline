# Dockerfile for ENCODE-DCC mirna-seq-pipeline
# base on ubuntu 20.04
FROM ubuntu@sha256:10cbddb6cf8568f56584ccb6c866203e68ab8e621bb87038e254f6f27f955bbe
MAINTAINER Otto Jolanki

RUN apt-get update && apt-get install -y \
    python3-dev \
    python3-pip \
    wget \
    libkrb5-3 \
    git \
    pigz \
    #samtools dependencies
    libbz2-dev \
    libncurses5-dev \
    libcurl4

RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install STAR/Samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install && rm ../zlib-1.2.11.tar.gz
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install && rm ../xz-5.2.3.tar.gz

RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz
RUN cd STAR-2.5.1b && make STAR && rm ../2.5.1b.tar.gz
ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

# Install wigToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && chmod +x wigToBigWig

RUN pip3 install cutadapt==3.4

RUN git clone --branch 1.9 --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch 1.9 --single-branch git://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

RUN pip3 install qc-utils==0.1.1
RUN pip3 install pandas==1.3.2 scipy==1.7.1

RUN mkdir -p mirna-seq-pipeline/src
COPY /src mirna-seq-pipeline/src
ENV PATH="/software/mirna-seq-pipeline/src:${PATH}"
