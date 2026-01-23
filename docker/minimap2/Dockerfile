FROM staphb/minimap2:2.28

# Install samtools
RUN apt-get update && apt-get install -y \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /data

