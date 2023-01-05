FROM ubuntu:xenial

# for easy upgrade later. ARG variables only persist during build time.
ARG K2VER="2.1.2"
ARG BVER="2.8"
ARG KRVER="2.8.1"

LABEL org.opencontainers.image.authors="tim.van.rillaer@hotmail.com"

RUN apt-get update && apt-get install --yes --no-install-recommends \
    wget \
    locales \
    nano \
    git \
    cmake \
    build-essential \
    default-jre \
    gcc-multilib \
    perl \
    g++ \
    rsync \
    curl \
    cpanminus \
    python \
    python-pip \
    r-base && \
    rm -rf /var/lib/apt/lists/* && apt-get autoclean

# perl modules
RUN cpanm Getopt::Std

# Install R libraries
RUN R -e 'install.packages("dplyr", repos="http://cloud.r-project.org/"); install.packages("randomForest",repos="http://cloud.r-project.org/")'

# Install Kraken2
RUN wget https://github.com/DerrickWood/kraken2/archive/v${K2VER}.tar.gz && \
 tar -xzf v${K2VER}.tar.gz && \
 rm -rf v${K2VER}.tar.gz && \
 cd kraken2-${K2VER} && \
 ./install_kraken2.sh . && \
 mkdir /data /kraken2-db

ENV PATH="$PATH:/kraken2-${K2VER}" \
    LC_ALL=C

# Install Bracken
RUN wget https://github.com/jenniferlu717/Bracken/archive/v${BVER}.tar.gz && \
 tar -xzf v${BVER}.tar.gz && \
 rm -rf v${BVER}.tar.gz && \
 cd Bracken-${BVER} && \
 chmod +x install_bracken.sh && \
 ./install_bracken.sh .

ENV PATH="$PATH:/Bracken-${BVER}" \
    LC_ALL=C

# Install Krona
RUN wget https://github.com/marbl/Krona/archive/v${KRVER}.tar.gz && \
tar -xzf v${KRVER}.tar.gz && \
 rm -rf v${KRVER}.tar.gz && \
 cd Krona-${KRVER}/KronaTools && \
 chmod +x install.pl && \
 ./install.pl . && \
 mkdir taxonomy && \
 ./updateTaxonomy.sh

# Install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
unzip fastqc_v0.11.9.zip && \
rm fastqc_v0.11.9.zip && \
cd FastQC && \
chmod +x fastqc && \
ln -s fastqc /usr/local/bin/fastqc

# Install multiqc
RUN pip install multiqc

