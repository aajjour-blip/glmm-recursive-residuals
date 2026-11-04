#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
command -v Rscript >/dev/null 2>&1 || { echo 'Rscript required'; exit 2; }
command -v python3 >/dev/null 2>&1 || { echo 'python3 required'; exit 2; }
mkdir -p results
# Parse first: never waste machine cycles on syntax errors.
Rscript -e 'parse(file="code/resolve_boundary.R"); cat("R_PARSE_PASS\n")' >results/parse.log
python3 code/preview_from_export.py
Rscript code/resolve_boundary.R 2>&1 | tee results/boundary_run.log
python3 - <<'PYCODE'
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib
root=Path('results');p=Path('GLMM_1725547_CBPP_BOUNDARY_RESOLUTION_RESULTS.zip')
files=sorted(x for x in root.rglob('*') if x.is_file())
with ZipFile(p,'w',ZIP_DEFLATED,compresslevel=9) as z:
 for file in files:z.write(file,str(file))
with ZipFile(p) as z:
 assert z.testzip() is None
 print('ZIP_INTEGRITY_PASS',len(z.namelist()))
print('RESULTS_ZIP',p.resolve())
print('RESULTS_SHA256',hashlib.sha256(p.read_bytes()).hexdigest())
PYCODE
