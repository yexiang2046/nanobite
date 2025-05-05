FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('optparse', 'tidyverse', 'pheatmap', 'RColorBrewer'), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('DESeq2', 'BiocParallel'), update=FALSE, ask=FALSE)"

# Set working directory
WORKDIR /data

# Create output directories
RUN mkdir -p /data/results/deseq2/results \
    /data/results/deseq2/plots

# Set environment variables
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Copy R scripts
COPY bin/prepare_deseq2.R /usr/local/bin/
COPY bin/run_deseq2.R /usr/local/bin/

# Make scripts executable
RUN chmod +x /usr/local/bin/prepare_deseq2.R \
    /usr/local/bin/run_deseq2.R

# Set default command
CMD ["R"] 