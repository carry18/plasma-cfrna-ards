"""Check all final V26 deliverables and statistical invariants; emit a reviewable report."""
import csv,json,hashlib,sys,re
from pathlib import Path
P=Path('results'); checks=[]
def check(name,ok,detail=''):
    checks.append(dict(name=name,passed=bool(ok),detail=detail))
def rows(p):
    with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
for n in range(1,7):
    for ext in ['svg','png','pdf','tiff']:
        f=P/f'figures/Figure{n}/Figure{n}.{ext}';check(str(f.relative_to(P)),f.exists() and f.stat().st_size>100)
for n in [1,2]:
    for ext in ['svg','png','pdf','tiff']:
        f=P/f'supplementary_figures/FigureS{n}.{ext}';check(str(f.relative_to(P)),f.exists() and f.stat().st_size>100)
expected=['TableS1_DE_summary','TableS2_Hallmark_all_contrasts','TableS3a_module_sizes','TableS3b_all_patient_overlap','TableS3c_module_effects','TableS3d_M02_resampling_summary','TableS3e_M02_resampling_draws','TableS4a_clinical_associations','TableS4b_platelet_adjustment','TableS5a_external_score','TableS5b_external_kME']
for name in expected:
    f=P/f'supplementary_tables/{name}.csv';check(name,f.exists() and bool(rows(f)))
md=rows(P/'analysis/included_samples_patients79.csv')
check('79 independent patients',len(md)==79 and len({r['Sample'] for r in md})==79)
check('groups A/B/C=29/17/33',[sum(r['group']==g for r in md) for g in 'ABC']==[29,17,33])
check('uploaded sample IDs', all(re.fullmatch(r'CFPL[0-9]{3}', r['Sample']) for r in md))
sizes=rows(P/'supplementary_tables/TableS3a_module_sizes.csv');check('M02 1159 genes',any(r['module']=='M02' and int(r['size'])==1159 for r in sizes))
clinical=rows(P/'supplementary_tables/TableS4a_clinical_associations.csv');check('30 clinical tests',len(clinical)==30)
platelet=next(r for r in clinical if r['signature']=='M02' and r['endpoint']=='platelet')
check('platelet n=40 and nominal P',int(platelet['n'])==40 and abs(float(platelet['P'])-.00412)<.00001)
adjustment=next(r for r in rows(P/'supplementary_tables/TableS4b_platelet_adjustment.csv') if r['program']=='M02_module')
check('platelet-adjusted M02 P=0.0674',abs(float(adjustment['adjusted_P'])-.0674)<.0001)
external=rows(P/'supplementary_tables/TableS5a_external_score.csv')[0]
check('external 684 genes, beta=0.816',int(external['measured_genes'])==684 and abs(float(external['estimate'])-.816)<.001)
draws=rows(P/'supplementary_tables/TableS3e_M02_resampling_draws.csv')
check('350 resampling draws',len(draws)==350 and [sum(r['mode']==mode and float(r['fraction'])==frac for r in draws) for mode,frac in [('fixed',.8),('fixed',.9),('reselected',.8),('loo',1)]]==[100,100,100,50])
check('all B/C directions retained',all(float(r['direction_retained'])==1 for r in draws))
t1=rows(P/'Table1.csv')
check('Table 1 severity and ICU rows distinct',any(r['Characteristic']=='ICU stay days, median (IQR)' for r in t1) and any(r['Characteristic']=='APACHE II score, median (IQR)' for r in t1))
# All implementation sources are included in a checksum inventory.
manifest=[]
for f in sorted(Path('scripts').rglob('*')):
    if f.is_file() and '__pycache__' not in f.parts:manifest.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(P/'code_checksums.json').write_text(json.dumps(manifest,indent=2)+'\n')
report={'status':'passed' if all(c['passed'] for c in checks) else 'failed','checks':checks,'scope':'Final manuscript statistical and artifact invariants; consult logs and validation_notes.md for recomputation/cache provenance.'}
(P/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(f"{report['status']}: {sum(c['passed'] for c in checks)}/{len(checks)} checks")
for c in checks:
    if not c['passed']:print(c['name'],c['detail'])
sys.exit(0 if report['status']=='passed' else 1)
