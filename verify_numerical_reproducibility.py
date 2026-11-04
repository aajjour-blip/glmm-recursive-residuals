#!/usr/bin/env python3
"""Evidence-based GLMM manuscript 1725547 numeric checks. Never modifies archived results.

Usage:
    python3 verify_numerical_reproducibility.py
    python3 verify_numerical_reproducibility.py --strict-exact

Default: independently reproduce archived-output checks and disclose observed
stochastic-run exceptions. Strict: nonzero when any new-run mismatch remains.
"""
from pathlib import Path
import csv, json, hashlib, subprocess, tempfile, shutil, sys, collections
ROOT=Path(__file__).resolve().parent
EVID=ROOT/'verification_evidence/generated'
report={'manuscript':'1725547','scope':'Archived intermediate-data reconstruction and Mac-run comparison; not exact rerun certification','checks':{}, 'exceptions':[]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok, msg):
 if not ok:raise AssertionError(msg)
def verify_file(a,b,mode='bytes'):
 require(a.exists(),f'Input absent: {a}')
 require(b.exists(),f'Reference absent: {b}')
 if mode=='bytes':return sha(a)==sha(b)
 if mode=='whitespace':return ''.join(a.read_text().split())==''.join(b.read_text().split())
 raise ValueError(mode)
# Independent regeneration from archival CSVs using the official Python table generator.
with tempfile.TemporaryDirectory(prefix='glmm_rr_tables_') as tmp:
 t=Path(tmp)
 (t/'data').symlink_to(ROOT/'data',target_is_directory=True)
 (t/'tables').mkdir()
 command=[sys.executable,str(ROOT/'code/07_format_tables.py')]
 outcome=subprocess.run(command,cwd=t,capture_output=True,text=True)
 require(outcome.returncode==0,f'Table generator failed: {outcome.stdout} {outcome.stderr}')
 x=sorted((t/'tables').glob('*.tex'))
 require(len(x)==5,f'Expected five generated tables, got {len(x)}')
 require(all(verify_file(f,ROOT/'tables'/f.name) for f in x),'Archived data-to-table reproduction mismatch')
 report['checks']['archived_csv_to_latex_tables']={'status':'PASS_EXACT','count':5,'names':[f.name for f in x]}
# Manuscript LaTeX from already completed Mac scripts 08 and 09.
exact=[]; format_only=[]
for id in ('08','09'):
 fs=sorted((EVID/id/'tables').glob('*.tex'))
 require(fs,f'No generated LaTeX for workflow {id}')
 for f in fs:
  r=ROOT/'tables'/f.name
  if verify_file(f,r):exact.append(f.name)
  elif verify_file(f,r,'whitespace'):format_only.append(f.name)
  else:raise AssertionError('Numerical or content discrepancy in generated table: '+f.name)
require(len(exact)+len(format_only)==5,'Unexpected number of reporting tables')
report['checks']['mac_reporting_tables']={'status':'PASS_NUMERIC_CONTENT','exact_bytes':exact,'whitespace_only':format_only}
# Ablation raw CSV outputs and raw .tex from completed Mac script 06.
csv_files=sorted((EVID/'06'/'results/ablation_epilepsy/data').glob('*.csv'))
tex_files=sorted((EVID/'06'/'results/ablation_epilepsy/tables').glob('*.tex'))
require(len(csv_files)==6 and len(tex_files)==3,'Ablation evidence count incomplete')
require(all(verify_file(x,ROOT/'data/ablation_epilepsy'/x.name) for x in csv_files),'Ablation CSV mismatch')
raw_exact=[x.name for x in tex_files if verify_file(x,ROOT/'tables'/x.name)]
raw_different=[x.name for x in tex_files if x.name not in raw_exact]
report['checks']['mac_ablation']={'status':'PASS_CSV_EXACT_RAW_LATEX_PRESENTATION_DIFFERS','csv_count':6,'raw_latex_count':3,'raw_latex_byte_exact':raw_exact,'raw_latex_different':raw_different}
report['exceptions'].append('Three raw Mac ablation LaTeX tables differ from final manuscript LaTeX in layout/captions and one small p-value printed as 0.000 instead of <0.001. Final manuscript tables independently regenerate byte-for-byte from archived CSVs.')
# Original Bonferroni DHARMa result re-evaluation (uses archived replicate records).
outcome=subprocess.run([sys.executable,str(ROOT/'verify_archived_results.py')],cwd=ROOT,capture_output=True,text=True)
require(outcome.returncode==0,'Archived DHARMa audit failed: '+outcome.stdout+outcome.stderr)
report['checks']['dharma_archived_800_replicates']={'status':'PASS_INTERNAL','records':800}
# Verify outcome classifications by independently reading archived / reproduced CSVs.
outcome=subprocess.run([sys.executable,str(ROOT/'verification_evidence/independent_numeric_audit/reproduce_comparison.py')],cwd=ROOT,capture_output=True,text=True)
require(outcome.returncode==0,'Regenerated CSV audit failed: '+outcome.stdout+outcome.stderr)
with (ROOT/'verification_evidence/independent_numeric_audit/CSV_COMPARISON_DETAIL.csv').open(newline='') as f:rows=list(csv.DictReader(f))
cnt=dict(collections.Counter(row['classification'] for row in rows))
require(len(rows)==82 and cnt.get('NUMERIC_DIFFERENCE')==40,'Comparison count changed: '+str(cnt))
report['checks']['mac_vs_archived_csvs']={'status':'EXCEPTIONS','comparisons':len(rows),'classes':cnt}
report['exceptions'].append('Forty of 82 compared CSV outputs have numerical differences; unchanged omnibus detection outcomes do not establish equal p-values.')
report['exceptions'].append('Historical and causal replay evidence is now provided separately and supports comparator-order dependence in targeted prefix tests; it does not establish exact full-run B=1000 empirical reproduction.')
report['exceptions'].append('Full simulation rerun and independent regeneration of every main/supplementary cell are not demonstrated.')
# Verify full completed workflow execution evidence.
w=json.loads((ROOT/'verification_evidence/verification_summary.json').read_text())
require(w.get('all_executed_passed') is True and not w.get('unexecuted') and len(w['workflow_results'])==6,'Mac execution evidence incomplete')
require(all(r['return_code']==0 for r in w['workflow_results']),'Failed Mac workflow')
report['checks']['completed_mac_workflows']={'status':'PASS_EXECUTION','count':6,'generated_files':sum(r['result_files'] for r in w['workflow_results'])}
# Figure B1 pixel check was independently recorded; do not claim fresh visual rendering.
pdf=json.loads((ROOT/'verification_evidence/PDF_VISUAL_COMPARISON.json').read_text())
require(pdf.get('pixels_identical_at_1.5x') is True,'Ablation figure visual audit evidence missing')
report['checks']['ablation_pdf_prior_visual_comparison']={'status':'PASS_ARCHIVED_AUDIT','note':'Prior 1.5x pixel comparison; not newly rendered here'}
report['verdict']='PASS_ARCHIVED_RECONSTRUCTION_WITH_RECORDED_EXCEPTIONS'
print(json.dumps(report,indent=2))
if '--strict-exact' in sys.argv:
 print('STRICT_EXACT_NOT_MET: 40 pre-correction empirical CSV comparisons differ; full corrected-source B=1000 CSV regeneration has not been demonstrated.',file=sys.stderr)
 sys.exit(2)
