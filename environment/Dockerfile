# hash:sha256:37dbbfb5965d562dcc6138889f2e18924d951d97804c5b0c46a2b4a3ac31c442
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3-pip=9.0.1-2.3~ubuntu1.18.04.2 \
        python3-setuptools=39.0.1-2 \
        python3-wheel=0.30.0-0.2 \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir \
    jupyterlab==2.2.4
