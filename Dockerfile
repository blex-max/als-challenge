# ── build stage ───────────────────────────────────────────────────────────────
# python:3.12 satisfies requires-python = ">=3.12" and is stable on Docker Hub.
FROM python:3.12-slim-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        ninja-build \
        build-essential \
        pkg-config \
        libhts-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# scikit-build-core drives cmake → ninja, producing cfextract.so.
# All Python dependencies declared in pyproject.toml are installed too.
RUN pip install --no-cache-dir .

# ── runtime stage ─────────────────────────────────────────────────────────────
FROM python:3.12-slim-bookworm

# libhts3: the shared library cfextract.so links against at runtime.
# Build tools are not needed in the final image.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhts3 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/lib/python3.12/site-packages \
                    /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin/cfanalysis \
                    /usr/local/bin/cfanalysis

ENTRYPOINT ["cfanalysis"]
CMD ["--help"]
