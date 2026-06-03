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

# 5. Install Python Environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml .
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# 6. Copy Code
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# 7. Run Migrations
RUN micromamba run -n base python manage.py migrate

EXPOSE 8000

# 9. Start the server and run migrations at runtime
CMD ["micromamba", "run", "-n", "base", "bash", "-c", \
     "mkdir -p /app/data/cache /tmp/amber_django_cache; \
     chmod 666 db.sqlite3 2>/dev/null || true; \
     chmod -R 777 /app/data /tmp/amber_django_cache 2>/dev/null || true; \
     python manage.py migrate; python manage.py runserver 0.0.0.0:8000"]
