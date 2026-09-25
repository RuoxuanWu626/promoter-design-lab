#!/bin/bash
# Launch the promoter design server on a Randi compute node via Slurm,
# then print the exact SSH tunnel command to reach it from a laptop.
#
#   ./run_server.sh                 # submit to the default partition
#   ./run_server.sh -p express -t 6:00:00
#   ./run_server.sh --local         # run here instead (login node; short tests only)
#   ./run_server.sh --stop          # cancel a running server job
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARTITION="tier1q"
TIME="12:00:00"
CPUS=8
MEM="16G"
PORT="${PDG_PORT:-8765}"
# Default to the env that has torch, so the real Puffin checkpoint loads.
# Falls back to plain python3, where only the mock adapters are available.
TORCH_PY=/gpfs/data/zhou-lab/rxwu/settings/miniforge3/envs/alphagenome/bin/python
PYTHON="${PDG_PYTHON:-$([ -x "$TORCH_PY" ] && echo "$TORCH_PY" || echo python3)}"
JOBNAME="promoter-lab"
LOCAL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--partition) PARTITION="$2"; shift 2 ;;
    -t|--time)      TIME="$2"; shift 2 ;;
    -c|--cpus)      CPUS="$2"; shift 2 ;;
    -m|--mem)       MEM="$2"; shift 2 ;;
    --port)         PORT="$2"; shift 2 ;;
    --python)       PYTHON="$2"; shift 2 ;;
    --local)        LOCAL=1; shift ;;
    --stop)         scancel -n "$JOBNAME" -u "$USER" && echo "cancelled any '$JOBNAME' jobs"; exit 0 ;;
    -h|--help)      sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$ROOT/data/runs" "$ROOT/logs"

if [[ $LOCAL -eq 1 ]]; then
  echo "running on $(hostname) (login node — fine for a quick look, not for a class)"
  exec "$PYTHON" "$ROOT/backend/app.py" --port "$PORT" --open-port-scan
fi

SBATCH="$ROOT/logs/${JOBNAME}.sbatch"
cat > "$SBATCH" <<EOF
#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --partition=$PARTITION
#SBATCH --time=$TIME
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --output=$ROOT/logs/%x-%j.out
#SBATCH --error=$ROOT/logs/%x-%j.out

echo "node: \$(hostname)"
echo "port: $PORT"
cd "$ROOT/backend"
# Puffin's 601-wide deconvolutions are the bulk of the cost; give them the
# cores the job actually reserved.
export OMP_NUM_THREADS=$CPUS MKL_NUM_THREADS=$CPUS
exec $PYTHON app.py --port $PORT --host 0.0.0.0 --open-port-scan
EOF

JOBID=$(sbatch --parsable "$SBATCH")
echo "submitted job $JOBID to $PARTITION (${TIME}, ${CPUS} cpu, ${MEM})"
printf "waiting for it to start"

for _ in $(seq 1 120); do
  STATE=$(squeue -j "$JOBID" -h -o %T 2>/dev/null || echo "")
  [[ "$STATE" == "RUNNING" ]] && break
  [[ -z "$STATE" ]] && { echo; echo "job left the queue — check $ROOT/logs/"; exit 1; }
  printf "."; sleep 2
done
echo

NODE=$(squeue -j "$JOBID" -h -o %N)
LOG="$ROOT/logs/${JOBNAME}-${JOBID}.out"
for _ in $(seq 1 40); do
  grep -q "serving" "$LOG" 2>/dev/null && break
  sleep 1
done
ACTUAL_PORT=$(grep -oE 'bound [0-9.]+:[0-9]+' "$LOG" 2>/dev/null | tail -1 | cut -d: -f2)
ACTUAL_PORT="${ACTUAL_PORT:-$PORT}"

cat <<EOF

  server running on $NODE:$ACTUAL_PORT   (job $JOBID)
  log: $LOG

  From your laptop:

    ssh -N -L ${ACTUAL_PORT}:${NODE}:${ACTUAL_PORT} Randi

  then open  http://localhost:${ACTUAL_PORT}/

  Stop it with:  scancel $JOBID     (or ./run_server.sh --stop)
EOF
