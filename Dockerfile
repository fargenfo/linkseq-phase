FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="LinkSeq Phase [WIP]" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

RUN apt-get update -yqq && \
    apt-get install -yqq \
    unzip curl tar

WORKDIR /tmp

ENV HTSLIB_URL=https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
ENV HAPCUT2_URL=https://github.com/vibansal/HapCUT2/archive/c6481d5fd0618dc3e82b2eb8c2b4835d9a4f6da7.zip

RUN curl -L $HTSLIB_URL | tar -xj
RUN wget --quiet $HAPCUT2_URL && unzip -qq c6481d*.zip

RUN apt-get update -yqq && \
    apt-get install -yqq \
    gcc \
    make \
    libc-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev


RUN cd /tmp/htslib* && \
    ./configure --prefix=/usr/local/ && \
    make && \
    make install

RUN cd /tmp/HapCUT2* && \
    make && \
    make install

# Add path to libraries.
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/local/lib
