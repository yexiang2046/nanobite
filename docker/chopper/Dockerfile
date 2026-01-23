FROM debian:bookworm-slim

# Install curl and bash
RUN apt-get update && apt-get install -y curl bash \
    && curl -L https://github.com/wdecoster/chopper/releases/download/v0.12.0b/chopper-linux -o /usr/local/bin/chopper \
    && chmod +x /usr/local/bin/chopper \
    && apt-get purge -y curl \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Use bash as the entrypoint to allow piping and redirection
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["chopper --help"]