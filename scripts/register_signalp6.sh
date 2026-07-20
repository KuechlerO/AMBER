#!/bin/sh
# Register a DTU SignalP fast tarball into the container's signalp6 conda env.
# Usage:
#   ./scripts/register_signalp6.sh /path/to/signalp-6.0i.fast.tar.gz [container]
set -eu
TARBALL="${1:?Usage: $0 /path/to/signalp-6.0*.fast.tar.gz [container_name]}"
CONTAINER="${2:-git-userfolder-olis_amber_django_app-1}"
BASENAME="$(basename "$TARBALL")"

docker cp "$TARBALL" "$CONTAINER:/tmp/$BASENAME"
docker exec -u mambauser "$CONTAINER" sh -c "
  ENV=/opt/conda/envs/signalp6
  export PATH=\"\$ENV/bin:\$PATH\"
  export CONDA_PREFIX=\"\$ENV\"
  ln -sfn /tmp/$BASENAME /tmp/signalp-6.0h.fast.tar.gz
  signalp6-register /tmp/signalp-6.0h.fast.tar.gz
  signalp6 -h | head -5
"
