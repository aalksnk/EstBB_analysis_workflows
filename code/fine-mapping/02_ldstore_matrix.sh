#!/bin/bash
#SBATCH --job-name=ldstore_regions
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --array=1-120             # set to number of z files
#SBATCH --output=logs/ldstore_%A_%a.out
#SBATCH --error=logs/ldstore_%A_%a.err

set -euo pipefail

#########################
# MODULES / CONTAINER
#########################
module load any/jdk/1.8.0_265
module load squashfs/4.4
module load any/singularity/3.5.3 2>/dev/null || true
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true

IMAGE="/path/quay.io-idarahu-ldmatrix:v0.1.img"

#########################
# PATH SETTINGS
#########################

# Directory with all .z files
ZDIR="/path"

# Directories for bgen/bgi/sample/inclusion
BGEN_DIR="/path"            # expects chr<CHR>.bgen
SAMPLE_DIR="/pathB"       # expects chr<CHR>.sample
GLOBAL_INCL="/path"      # expects one per region, e.g. chr1_.._...(.incl/.tsv/.txt)

# Work directory 
WORKDIR="/path"
mkdir -p "$WORKDIR"


LD_OUTDIR="/path"
BCOR_OUTDIR="/path"
mkdir -p "$LD_OUTDIR" "$BCOR_OUTDIR"

if [[ ! -f "$GLOBAL_INCL" ]]; then
  echo "ERROR: Global inclusion file not found: $GLOBAL_INCL"
  exit 1
fi

set -euo pipefail

# Where to collect failures (choose a shared path)
FAIL_TABLE="${LD_OUTDIR}/ldstore_failures.tsv"
LOCKFILE="${FAIL_TABLE}.lock"

init_fail_table() {
  mkdir -p "$(dirname "$FAIL_TABLE")"
  # Ensure header exactly once under lock
  {
    flock -x 200
    if [[ ! -f "$FAIL_TABLE" ]]; then
      printf "job_id\tarray_task_id\tidx\tzfile\tzbase\tregion\tchr\tstart\tend\texit_code\treason\twhen\n" \
        > "$FAIL_TABLE"
    fi
  } 200>"$LOCKFILE"
}

record_failure() {
  local exit_code="$1"
  local reason="$2"
  local when
  when="$(date -Is)"

  # Use safe expansions in case failure happens early
  local job_id="${SLURM_JOB_ID:-manual}"
  local array_id="${SLURM_ARRAY_TASK_ID:-NA}"
  local idx="${IDX:-NA}"
  local zfile="${ZFILE:-NA}"
  local zbase="${ZBASE:-NA}"
  local region="${REGION:-NA}"
  local chr="${CHR:-NA}"
  local start="${START:-NA}"
  local end="${END:-NA}"

  {
    flock -x 200
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$job_id" "$array_id" "$idx" "$zfile" "$zbase" "$region" \
      "$chr" "$start" "$end" "$exit_code" "$reason" "$when" \
      >> "$FAIL_TABLE"
  } 200>"$LOCKFILE"
}

init_fail_table

# Log any non-zero exit automatically
trap '
  ec=$?
  if [[ $ec -ne 0 ]]; then
    reason="error"
    # Common-ish OOM/kill signals as seen by the shell
    if [[ $ec -eq 137 || $ec -eq 134 || $ec -eq 9 ]]; then
      reason="OOM suspected / killed"
    fi
    record_failure "$ec" "$reason"
  fi
  exit $ec
' EXIT



#########################
# PICK THIS JOB'S Z FILE
#########################

mapfile -t ZFILES < <(ls "${ZDIR}"/*.z)

IDX=$((SLURM_ARRAY_TASK_ID - 1)) || IDX=0

if (( IDX < 0 || IDX >= ${#ZFILES[@]} )); then
  echo "Array index ${SLURM_ARRAY_TASK_ID} out of range (have ${#ZFILES[@]} z-files)"
  exit 1
fi

ZFILE="${ZFILES[$IDX]}"
ZBASE="$(basename "$ZFILE" .z)"

echo "Processing Z file: $ZFILE"

#########################
# DERIVE CHR + REGION
#########################

# Expect EXACT base name:
# reference_LD_<CHR>_<LDSTART>_<LDEND>
# e.g. reference_LD_3_31103250_32103250

if [[ "$ZBASE" =~ ^reference_LD_([0-9]+)_([0-9]+)_([0-9]+)$ ]]; then
    RAWCHR="${BASH_REMATCH[1]}"
    START="${BASH_REMATCH[2]}"
    END="${BASH_REMATCH[3]}"

    CHR="$RAWCHR"
    REGION="${CHR}:${START}-${END}"
else
    echo "ERROR: Expected ZBASE format 'reference_LD_<chr>_<ldstart>_<ldend>' but got: $ZBASE"
    exit 1
fi

echo "Parsed CHR=${CHR}, REGION=${REGION}"


#########################
# SKIP IF OUTPUTS EXIST
#########################

OUT_BCOR="${BCOR_OUTDIR}/${REGION}.bcor"
OUT_LD="${LD_OUTDIR}/${REGION}.ld.gz"

if [[ -s "$OUT_BCOR" && -s "$OUT_LD" ]]; then
  echo "Outputs already exist for ${REGION}:"
  echo "  $OUT_BCOR"
  echo "  $OUT_LD"
  echo "Skipping."
  exit 0
fi


#########################
# FIND MATCHING FILES
#########################

# BGEN/BGI for that chromosome
BGEN="${BGEN_DIR}_chr${CHR}.bgen"
BGI="${BGEN}.bgi"

if [[ ! -f "$BGEN" ]]; then
    echo "ERROR: BGEN not found: $BGEN"
    exit 1
fi

if [[ ! -f "$BGI" ]]; then
    echo "ERROR: BGI not found: $BGI"
    exit 1
fi

# SAMPLE for that chromosome
SAMPLE_FILE="${SAMPLE_DIR}_chr${CHR}.sample"
if [[ ! -f "$SAMPLE_FILE" ]]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi



#########################
# PREP PER-JOB WORKDIR
#########################

JOBDIR="${WORKDIR}/job_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$JOBDIR"

# copy inputs into job-local dir
cp -f "$ZFILE"       "${JOBDIR}/${REGION}.z"
cp -f "$BGEN"        "${JOBDIR}/"
cp -f "$BGI"         "${JOBDIR}/"
cp -f "$SAMPLE_FILE" "${JOBDIR}/"

# Copy global inclusion file as incl.incl
cp -f "$GLOBAL_INCL" "${JOBDIR}/incl.incl"

# sample count from inclusion list
NSAMPLES=$(wc -l < "${JOBDIR}/incl.incl")

BGEN_BASE="$(basename "$BGEN")"
BGI_BASE="$(basename "$BGI")"
SAMPLE_BASE="$(basename "$SAMPLE_FILE")"

#########################
# MASTER FILE
#########################

cat > "${JOBDIR}/master.txt" <<EOF
z;bgen;bgi;bcor;ld;n_samples;sample;incl
/work/${REGION}.z;/work/${BGEN_BASE};/work/${BGI_BASE};/work/${REGION}.bcor;/work/${REGION}.ld;${NSAMPLES};/work/${SAMPLE_BASE};/work/incl.incl
EOF

#########################
# RUN LDSTORE
#########################

THREADS=${SLURM_CPUS_PER_TASK:-1}

singularity exec -B "${JOBDIR}:/work" --pwd /work "$IMAGE" bash -lc "
  set -euo pipefail

  echo '--- working directory ---'
  pwd

  echo '--- master.txt ---'
  cat master.txt

  echo '--- removing old possibly failed outputs ---'
  rm -f '${REGION}.bcor' '${REGION}.ld' '${REGION}.ld.gz'

  echo '--- creating .bcor ---'
  ldstore \
    --in-files master.txt \
    --write-bcor \
    --read-only-bgen \
    --n-threads ${THREADS}

  echo '--- checking newly created .bcor ---'
  if [ ! -s '${REGION}.bcor' ]; then
    echo 'ERROR: ${REGION}.bcor was not created or is empty'
    ls -lh
    exit 1
  fi

  ls -lh '${REGION}.bcor'

  echo '--- converting .bcor to text LD ---'
  ldstore \
    --in-files master.txt \
    --bcor-to-text

  echo '--- checking text LD output ---'
  if [ ! -s '${REGION}.ld' ]; then
    echo 'ERROR: ${REGION}.ld was not created or is empty'
    ls -lh
    exit 1
  fi

  ls -lh '${REGION}.ld'
"

# Compress LD
gzip -f "${JOBDIR}/${REGION}.ld"
#########################
# SAVE OUTPUTS
#########################

# Use REGION-based filenames; adjust to ZBASE if you prefer
cp -f "${JOBDIR}/${REGION}.bcor"    "${BCOR_OUTDIR}/${REGION}.bcor"
cp -f "${JOBDIR}/${REGION}.ld.gz"   "${LD_OUTDIR}/${REGION}.ld.gz"

echo "Done: ${ZBASE} -> ${BCOR_OUTDIR}/${REGION}.bcor and ${LD_OUTDIR}/${REGION}.ld.gz"
