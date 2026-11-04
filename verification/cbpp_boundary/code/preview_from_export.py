#!/usr/bin/env python3
"""Independent CSV-precision preview only; does not replace Mac full-precision validation."""
import csv,math,collections
from pathlib import Path
here=Path(__file__).resolve().parents[1]
r=list(csv.DictReader((here/'reference/singular_refit_row_predictions.csv').open()))
fit={int(x['replicate']):float(x['sigma']) for x in csv.DictReader((here/'reference/singular_refit_fit_and_prediction_diagnostics.csv').open())}
by=collections.defaultdict(list)
for x in r:by[int(x['replicate'])].append(x)
with (here/'results/csv_precision_preview.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['replicate','sigma','preview_Tinf','preview_T2','note'])
 for i,rows in sorted(by.items()):
  if fit[i]==0:
   grp=collections.defaultdict(list)
   for x in rows:grp[x['period_num']].append(float(x['pred_prob']))
   means={k:sum(v)/len(v) for k,v in grp.items()}
   keys=[means[x['period_num']] for x in rows]
  else:keys=[float(x['pred_prob']) for x in rows]
  rr=[float(x['rr']) for x in rows]
  idx=sorted(range(len(rr)),key=lambda j:(keys[j],j))
  c=0.;vals=[]
  for j in idx:c+=rr[j];vals.append(c/math.sqrt(len(rr)))
  w.writerow([i,fit[i],max(map(abs,vals)),sum(x*x for x in vals)/len(vals),'CSV-precision approximation only'])
 print('INDEPENDENT_PREVIEW_PASS',len(by),'singular refits')
