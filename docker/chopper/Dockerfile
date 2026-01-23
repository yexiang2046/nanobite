FROM debian:bookworm-slim

# Install curl to download the binary
RUN apt-get update && apt-get install -y curl \
    && curl -L https://github.com/wdecoster/chopper/releases/download/v0.12.0b/chopper-linux -o /usr/local/bin/chopper \
    && chmod +x /usr/local/bin/chopper \
    && apt-get purge -y curl \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Set the entrypoint to run chopper
ENTRYPOINT ["chopper"]
CMD ["--help"]