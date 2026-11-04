#!/usr/bin/env python3
"""Read-only numerical comparison of archived and regenerated outputs."""
from pathlib import Path
import csv,json,math,hashlib
root=Path(__file__).resolve().parent
report=json.loads((root/'verification/empirical_comparison_summary.json').read_text())
assert report['RR_files']==69 and report['DHARMa_files']==7
rr=dh=rrdif=dhdif=changes=0
for entry in report['details']:
    a=root/entry['archived']; b=root/entry['regenerated']
    if not a.is_file() or not b.is_file():raise SystemExit('MISSING_INPUT '+str(a)+' | '+str(b))
    with a.open(newline='',encoding='utf-8-sig') as f:x=list(csv.reader(f))
    with b.open(newline='',encoding='utf-8-sig') as f:y=list(csv.reader(f))
    assert len(x)==len(y) and all(len(m)==len(n) for m,n in zip(x,y)),entry['file']
    count=0
    for i,(m,n) in enumerate(zip(x,y)):
        for j,(s,t) in enumerate(zip(m,n)):
            if s==t:continue
            try:
                fs,ft=float(s),float(t)
                equal=math.isclose(fs,ft,rel_tol=0,abs_tol=1e-7)
                if not equal and entry['kind']=='DHARMa' and 'p' in x[0][j].lower():
                    changes+=((fs<.05)!=(ft<.05))
            except (TypeError,ValueError):equal=False
            if not equal:count+=1
    assert count==entry['different_cells'],(entry['file'],count,entry['different_cells'])
    if entry['kind']=='RR_OR_OTHER':rr+=1;rrdif+=count
    else:dh+=1;dhdif+=count
print('FINAL_RR',rr,'files',rrdif,'numeric/text discrepancies')
print('FINAL_DHARMA',dh,'files',dhdif,'differences',changes,'5% decision changes')
assert (rr,rrdif,dh,dhdif,changes)==(69,0,7,3,0)
print('FINAL_VERIFY_PASS')
