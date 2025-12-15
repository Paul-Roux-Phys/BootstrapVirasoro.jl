#!/usr/bin/env bash

# ---- configuration ----
REMOTE_USER="roux"
REMOTE_HOST="thanos"
REMOTE="${REMOTE_USER}@${REMOTE_HOST}"

# ---- argument parsing ----
CPUS="${1:-4}"

if ! [[ "$CPUS" =~ ^[0-9]+$ ]] || [[ "$CPUS" -lt 1 ]]; then
    echo "❌ Invalid CPU count: $CPUS" >&2
    exit 1
fi

echo "🔌 Requesting SLURM node on $REMOTE_HOST ($CPUS CPUs)..." >&2

REMOTE_COMMAND=$(cat <<EOF
set -e

OUTFILE=\$(mktemp)

(
  srun --partition=C-Short --account=thanos \
       -N 1 -n 1 --cpus-per-task=${CPUS} --mem=32G --time=12:00:00 \
       bash -lc '
         echo "NODE=\$(hostname)"
         /home/roux/.julia/conda/3/x86_64/bin/jupyter notebook --no-browser --ip=0.0.0.0 --port=8888
       '
) &> "\$OUTFILE" &

for i in {1..30}; do
    NODE=\$(grep "^NODE=" "\$OUTFILE" | sed 's/NODE=//')
    TOKEN=\$(grep -o "token=[a-f0-9]\+" "\$OUTFILE" | head -n1 | sed 's/token=//')
    if [[ -n "\$NODE" && -n "\$TOKEN" ]]; then
        echo "\$NODE" > "\$OUTFILE.final"
        echo "\$TOKEN" >> "\$OUTFILE.final"
        break
    fi
    sleep 1
done

cat "\$OUTFILE.final" 2>/dev/null || echo "__ERROR__"
EOF
)

RESULT=$(ssh "$REMOTE" "$REMOTE_COMMAND")

if [[ "$RESULT" == "__ERROR__" ]]; then
    echo "❌ Failed to extract node/token remotely." >&2
    exit 1
fi

NODE=$(echo "$RESULT" | sed -n '1p')
TOKEN=$(echo "$RESULT" | sed -n '2p')

if [[ -z "$NODE" || -z "$TOKEN" ]]; then
    echo "❌ Bad result:" >&2
    echo "$RESULT" >&2
    exit 1
fi

echo "🖥️  Compute node: $NODE"
echo "🔑  Token: $TOKEN"

LOCAL_PORT=$(python3 -c "import socket; s=socket.socket(); s.bind(('',0)); print(s.getsockname()[1]); s.close()")

echo "🧵 Opening SSH tunnel (${LOCAL_PORT} → $NODE:8888)…"
autossh -M 0 -f -N -o ControlMaster=no -L ${LOCAL_PORT}:${NODE}:8888 "$REMOTE"

URL="http://localhost:${LOCAL_PORT}/?token=${TOKEN}"
echo -n "$URL" | pbcopy

echo "📋 Copied to clipboard:"
echo "    $URL"
