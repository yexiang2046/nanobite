FROM nvidia/cuda:12.1.0-runtime-ubuntu22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install required dependencies
RUN apt-get update && apt-get install -y \
    wget \
    tar \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Create directories
WORKDIR /opt/dorado

# Download and install Dorado v0.9.5
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.5-linux-x64.tar.gz \
    && tar -xzvf dorado-0.9.5-linux-x64.tar.gz \
    && rm dorado-0.9.5-linux-x64.tar.gz

# Add Dorado to PATH
ENV PATH="/opt/dorado/dorado-0.9.5-linux-x64/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/dorado/dorado-0.9.5-linux-x64/lib:${LD_LIBRARY_PATH}"

# Set working directory for the application
WORKDIR /workspace

# Default command
CMD ["bash"] 