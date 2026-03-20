ARG QIIME_TAG=2026.1
FROM quay.io/qiime2/tiny:${QIIME_TAG}

# Build from the Plugins workspace root:
# docker build -f q2-bowtie2/Dockerfile -t q2-bowtie2 .
COPY q2-types-variant /plugins/q2-types-variant
COPY q2-bowtie2 /plugins/q2-bowtie2

RUN conda install -y -c conda-forge -c bioconda bowtie2 samtools && \
    conda clean -afy && \
    python -m pip install --no-cache-dir hatchling && \
    python -m pip install --no-cache-dir --no-build-isolation --no-deps /plugins/q2-types-variant && \
    python -m pip install --no-cache-dir --no-build-isolation --no-deps /plugins/q2-bowtie2 && \
    qiime dev refresh-cache
