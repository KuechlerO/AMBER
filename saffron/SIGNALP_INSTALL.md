SignalP 6.0 sidecar environment (SAFFRON)

Why a separate env?
  Django runs on Python 3.12 (see environment.yml). SignalP 6.0 needs
  Python ≤ 3.10, so it lives in its own micromamba env.

1) Create the conda env (stub)
  micromamba create -y -n signalp6 -f environment.signalp.yml

  This installs `predector::signalp6` (latest channel build is **6.0h**).
  That package is a **license stub** until you register a DTU tarball.

2) Download the licensed fast package from DTU
  https://services.healthtech.dtu.dk/services/SignalP-6.0h/9-Downloads.php#
  (or the matching SignalP-6.0i download page if that is what DTU issued you)

  Preferred filename for the stub: signalp-6.0h.fast.tar.gz
  A 6.0i fast tarball also works with the 6.0h stub if you symlink the name
  (the archive still contains signalp6_fast/ as expected).

3) Register inside the Docker container (important details!)
  The register script needs CONDA_PREFIX pointed at the **signalp6** env
  (not the base env). Without that it looks for unregister.sh in the wrong place.

  Example with a file already copied to /tmp:

    docker cp signalp-6.0i.fast.tar.gz git-userfolder-olis_amber_django_app-1:/tmp/

    docker exec -u mambauser git-userfolder-olis_amber_django_app-1 sh -c '
      ENV=/opt/conda/envs/signalp6
      export PATH="$ENV/bin:$PATH"
      export CONDA_PREFIX="$ENV"
      ln -sfn /tmp/signalp-6.0i.fast.tar.gz /tmp/signalp-6.0h.fast.tar.gz
      signalp6-register /tmp/signalp-6.0h.fast.tar.gz
      signalp6 -h
    '

  Success looks like: “signalp6 is now fully installed!” and `signalp6 -h`
  prints a usage message (not the license stub text).

  Note: plain `signalp6` is NOT on PATH in the base env. Always use
    /opt/conda/envs/signalp6/bin/signalp6
  or `export PATH=/opt/conda/envs/signalp6/bin:$PATH`.
  Django uses SIGNALP6_BIN=/opt/conda/envs/signalp6/bin/signalp6.

4) Point Django at the CLI
  SIGNALP6_BIN=/opt/conda/envs/signalp6/bin/signalp6
  SIGNALP_CACHE_DIR=/tmp/saffron_signalp

  Until step 3 succeeds, SAFFRON detects the stub and uses the mock predictor.

Force mock (tests / CI)
  SIGNALP_MOCK=1

Optional: auto-register on container start
  Mount the tarball and set SIGNALP6_TARBALL to its path (see Dockerfile CMD).
