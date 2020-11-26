FROM rocker/r-ver:3.5.3


RUN apt-get update --fix-missing -qq && \
    apt-get install -y -q \
    vim \
    git \
    python3 \
    python3-pip \
    python \
    python-pip \
    libz-dev \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libpng-dev \
    libjpeg-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libgl-dev \
    libgsl-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
#
RUN R -e 'install.packages(c("BiocManager", "devtools", "argparse"))'
#
RUN mkdir -p /usr/local/lib/R/site-library/
ADD latticeExtra /usr/local/lib/R/site-library/latticeExtra
#
RUN R -e 'install.packages("XML", repos = "http://www.omegahat.net/R")'
#
RUN R -e 'install.packages("Hmisc")'
#
RUN R -e 'BiocManager::install(c("rhdf5", "EnrichmentBrowser", "DRIMSeq", "tximport"))'
RUN R -e 'BiocManager::install("DEXSeq")'
#
RUN pip install numpy
#
ADD software /home/software
ADD data /home/data
ADD scripts /home/scripts
#
RUN ln /usr/lib/x86_64-linux-gnu/libgsl.so.19 /usr/lib/x86_64-linux-gnu/libgsl.so.0

