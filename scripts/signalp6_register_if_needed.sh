#!/bin/sh
# Register DTU SignalP 6.0 inside the signalp6 conda env when a tarball is mounted.
# Called on container start (see Dockerfile CMD / docker-compose SIGNALP6_TARBALL).
set -eu

ENV=/opt/conda/envs/signalp6
export PATH="$ENV/bin:$PATH"
export CONDA_PREFIX="$ENV"

TARBALL="${SIGNALP6_TARBALL:-}"
if [ -z "$TARBALL" ] || [ ! -f "$TARBALL" ]; then
    exit 0
fi

if ! command -v signalp6-register >/dev/null 2>&1; then
    echo "WARN: signalp6-register not found; SAFFRON will use mock predictor." >&2
    exit 0
fi

# Skip if the CLI is already registered (fast restarts within the same container).
if out="$(signalp6 -h 2>&1)" || true; then
    case "$out" in
        *"has not been installed yet"*|*"signalp6-register"*) ;;
        *)
            echo "SignalP 6.0 already registered."
            exit 0
            ;;
    esac
fi

REGISTER_TARBALL="$TARBALL"
case "$(basename "$TARBALL")" in
    signalp-6.0i.fast.tar.gz)
        ln -sfn "$TARBALL" /tmp/signalp-6.0h.fast.tar.gz
        REGISTER_TARBALL=/tmp/signalp-6.0h.fast.tar.gz
        ;;
esac

echo "Registering SignalP 6.0 from $REGISTER_TARBALL ..."
signalp6-register "$REGISTER_TARBALL"
echo "SignalP 6.0 registration complete."
