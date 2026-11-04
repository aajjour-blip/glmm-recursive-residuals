#!/usr/bin/env bash
set -Eeuo pipefail
cd "$(dirname "$0")"
DIR=fullscale_replay/runs/03/results/RR_guided_model_refinement_B1000_submission_results/data
[[ -d "$DIR" ]] || { echo '03 output directory is missing' >&2; exit 3; }
mkdir -p fullscale_replay/logs
python3 - <<'PY'
import csv
from pathlib import Path
p=Path('fullscale_replay/runs/03/results/RR_guided_model_refinement_B1000_submission_results/data')
for nm in ('Salamanders_NB2_env','Owls_ZINB_arrival','epil_NB2_refined'):
 f=p/(nm+'_RR_B1000_simulated_stats.csv')
 assert f.is_file(),str(f)
 with f.open(newline='') as h:n=sum(1 for _ in csv.DictReader(h))
 assert n==1000,(nm,n)
 print('RR_1000_PRESERVED',nm)
PY
export R_LIBS_USER="$PWD/fullscale_replay/historical_library${R_LIBS_USER:+:$R_LIBS_USER}"
Rscript -e 'stopifnot(as.character(packageVersion("DHARMa"))=="0.4.7");parse(file="recover_03_dharma.R");cat("RECOVERY_PARSE_OK\n")'
Rscript recover_03_dharma.R 2>&1 | tee fullscale_replay/logs/03_recovery.log
python3 - <<'PY'
import csv,sys
from pathlib import Path
sys.path.insert(0,'fullscale_replay')
from verify_fullscale import compare
base=Path('verification_evidence/independent_numeric_audit/CSV_COMPARISON_DETAIL.csv')
with base.open(newline='') as h:names=[r['file'] for r in csv.DictReader(h) if r['workflow']=='03']
assert len(names)==19,len(names)
new=Path('fullscale_replay/runs/03/results/RR_guided_model_refinement_B1000_submission_results/data')
old=Path('data/refinement');bad=[]
for name in names:
 a,b=new/name,old/name
 if not a.is_file() or not b.is_file():bad.append((name,'MISSING'));continue
 x=compare(a,b,1e-7)
 print('RECOVERY_COMPARE',name,x['status'],x.get('numeric_differences',0),x.get('text_differences',0))
 if x['status'] not in ('EXACT_BYTES','NUMERIC_AND_TEXT_MATCH'):bad.append((name,x['status']))
if bad:
 print('RECOVERY_03_UNRESOLVED',bad)
 sys.exit(2)
print('RECOVERY_03_19_OF_19_MATCH')
PY
