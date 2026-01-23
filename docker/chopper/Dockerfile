FROM rust:1.81-bookworm AS builder

# Install system dependencies required for building Rust crates with native dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libssl-dev \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone the repository and build chopper
WORKDIR /app
RUN git clone https://github.com/wdecoster/chopper.git .
RUN cargo build --release

# Final stage: Use a slim image for runtime
FROM debian:bookworm-slim

# Copy the built binary from the builder stage
COPY --from=builder /app/target/release/chopper /usr/local/bin/chopper

# Set the entrypoint to run chopper
ENTRYPOINT ["chopper"]
CMD ["--help"]