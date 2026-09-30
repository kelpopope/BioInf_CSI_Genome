FROM python:3.12.14-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /project
COPY requirements.lock.txt /tmp/requirements.lock.txt
RUN pip install --no-cache-dir -r /tmp/requirements.lock.txt
CMD ["snakemake", "--cores", "2"]
