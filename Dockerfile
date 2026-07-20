# 1. Base Image
FROM mambaorg/micromamba:latest

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV HTTP_PROXY=http://proxy.charite.de:8080
ENV HTTPS_PROXY=http://proxy.charite.de:8080

# 2. System installs (as root)
USER root
RUN apt-get update && apt-get install -y libpq-dev gcc && rm -rf /var/lib/apt/lists/*

# 3. Setup App Directory and PERMISSIONS (as root)
WORKDIR /app

# We create the data folder AND ensure mambauser owns /app so pip can write temp files
RUN chown -R $MAMBA_USER:$MAMBA_USER /app

# 4. Setup App
USER $MAMBA_USER

# 5. Install Python Environment (Django app)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml .
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# 5b. Optional SignalP 6.0 sidecar (Python 3.10) — see environment.signalp.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.signalp.yml .
RUN micromamba create -y -n signalp6 -f environment.signalp.yml && \
    micromamba clean --all --yes || \
    echo "WARN: SignalP env install failed; SAFFRON will use mock predictor"

ENV SIGNALP6_BIN=/opt/conda/envs/signalp6/bin/signalp6
ENV SIGNALP_CACHE_DIR=/tmp/saffron_signalp

# 6. Copy Code
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# 7. Run Migrations
RUN micromamba run -n base python manage.py migrate

EXPOSE 8000

# SignalP 6.0: second micromamba env `signalp6` from environment.signalp.yml.
# SIGNALP6_BIN defaults to /opt/conda/envs/signalp6/bin/signalp6.
# Mount signalp-6.0h.fast.tar.gz and set SIGNALP6_TARBALL to auto-register on start.
# See saffron/SIGNALP_INSTALL.md. Stub-only installs fall back to the mock predictor.

# 9. Start the server and run migrations at runtime
CMD ["micromamba", "run", "-n", "base", "bash", "-c", \
     "mkdir -p /app/data/cache /tmp/amber_django_cache /tmp/saffron_signalp; \
     chmod 666 db.sqlite3 2>/dev/null || true; \
     chmod -R 777 /app/data /tmp/amber_django_cache /tmp/saffron_signalp 2>/dev/null || true; \
     if [ -n \"$${SIGNALP6_TARBALL:-}\" ] && [ -f \"$${SIGNALP6_TARBALL}\" ]; then \
       /opt/conda/envs/signalp6/bin/signalp6-register \"$${SIGNALP6_TARBALL}\" || true; \
     fi; \
     python manage.py migrate; python manage.py runserver 0.0.0.0:8000"]
