#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PYTHON_BIN=${PYTHON:-python3}

if [[ -n "${SAGE:-}" ]]; then
  SAGE_BIN=$SAGE
elif command -v sage >/dev/null 2>&1; then
  SAGE_BIN=$(command -v sage)
elif [[ -x "/Applications/SageMath-10-9.app/Contents/Frameworks/Sage.framework/Versions/Current/local/bin/sage" ]]; then
  SAGE_BIN=/Applications/SageMath-10-9.app/Contents/Frameworks/Sage.framework/Versions/Current/local/bin/sage
else
  echo "SageMath 10.9 not found. Set SAGE=/path/to/sage." >&2
  exit 1
fi

"$SAGE_BIN" --version | grep -F "SageMath version 10.9" >/dev/null || {
  echo "The verification suite requires SageMath 10.9." >&2
  exit 1
}

"$SAGE_BIN" -python - <<'PY'
import sage.all
import snappy

if snappy.__version__ != "3.3.2":
    raise SystemExit(f"SnapPy 3.3.2 required, found {snappy.__version__}")
print("[PASS] shared SageMath 10.9 / SnapPy 3.3.2 runtime")
PY

cd "$SCRIPT_DIR"

"$PYTHON_BIN" riley_character.py
"$PYTHON_BIN" ptolemy_independent_check.py
"$PYTHON_BIN" riley_ptolemy_bridge.py
"$PYTHON_BIN" uniqueness_gcd.py
"$PYTHON_BIN" neumann_rogers_identity_check.py
"$PYTHON_BIN" dehn_wedge_reconstruction.py
"$PYTHON_BIN" anchor_integral_check.py

SAGE=$SAGE_BIN "$SAGE_BIN" -python byhand_flattening.py
SAGE=$SAGE_BIN "$SAGE_BIN" -python exactness_closure.py

echo "[PASS] Problem 3.1 independent verification suite"
