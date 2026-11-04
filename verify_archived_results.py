#!/usr/bin/env python3
import csv, math, pathlib, sys
p=pathlib.Path(__file__).parent/'data/completed_simulations/DHARMa_head_to_head/dharma_head_to_head_replicate_results.csv'
counts={}
with p.open(newline='') as f:
 for row in csv.DictReader(f):
  key=(row['scenario'],row['rep_id']); assert key not in counts, f'duplicate {key}'
  tests=[float(row[x]) for x in ('dharma_uniformity_p','dharma_dispersion_p','dharma_outlier_p')]
  assert all(math.isfinite(x) and 0<=x<=1 for x in tests)
  assert row['dharma_temporal_p'].strip().lower() in ('','na','nan'), 'temporal output unexpectedly populated'
  mi=min(tests); bo=min(1.,3*mi)
  assert abs(float(row['dharma_min_raw_p'])-mi)<1e-12
  assert abs(float(row['dharma_min_bonf_p'])-bo)<1e-12
  def boolean(v):return v.lower() in ('true','1','t')
  assert boolean(row['dharma_any_raw_detect'])==(mi<.05)
  assert boolean(row['dharma_any_bonf_detect'])==(bo<.05)
  counts[key]=row
assert len(counts)==800
expect={'null':(.07,.02,0.),'omitted_x2_gamma05':(.45,.03,0.),'omitted_x2_gamma10':(.99,.025,.005),'serial_dependence_rho10':(.375,.03,0.)}
for scenario, rates in expect.items():
 rows=[v for (s,r),v in counts.items() if s==scenario]; assert len(rows)==200
 rr=sum(boolean(v['rr_detect']) for v in rows)/200
 raw=sum(boolean(v['dharma_any_raw_detect']) for v in rows)/200
 bonf=sum(boolean(v['dharma_any_bonf_detect']) for v in rows)/200
 assert all(abs(a-b)<1e-12 for a,b in zip((rr,raw,bonf),rates))
 print(scenario,'n=200','RR=',rr,'DHARMa raw=',raw,'DHARMa Bonf=',bonf)
print('PASS: 800 records, no duplicate cells, 3-test p-values, combined rules and manuscript rates')
