FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="LinkSeq Phase [WIP]" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

RUN apt-get update -yqq && \
    apt-get install -yqq \
    unzip

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/linkseq-phase/bin:$PATH

# Pull the Nextflow pipeline.
# The following command makes sure the nextflow pul command isn't cached.
RUN CACHEBUST=$(date +%s)
RUN nextflow pull olavurmortensen/linkseq-phase
