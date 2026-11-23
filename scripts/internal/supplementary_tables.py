"""Generate Table 1 and supplementary tables from the analysis inputs and results."""
import csv,re,shutil,math,json
from pathlib import Path
ROOT=Path('results'); ANALYSIS=ROOT/'analysis'; TABLES=ROOT/'supplementary_tables'
TABLES.mkdir(parents=True,exist_ok=True)
def read(path):
    with Path(path).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def rows(sheet):
    data=list(sheet.values);return {i+2:dict(zip(data[0],r)) for i,r in enumerate(data[1:])}
def blank(v):return v is None or str(v).strip() in {'','NA','N/A','None'}
def num(v):
    m=re.search(r'(?<![\d-])-?\d+(?:\.\d+)?',str(v or ''));return float(m[0]) if m else None
def severity(v):
    m=re.search(r'APACHE\s*II\s*=\s*(\d+(?:\.\d+)?)',str(v or ''),re.I)
    if m:return float(m[1])
    ns=[float(x) for x in re.findall(r'(?<![\d-])\d+(?:\.\d+)?',str(v or '')) if not 1900<=float(x)<=2100]
    return ns[-1] if ns else None
def quantile(v,p):
    v=sorted(v);pos=(len(v)-1)*p;lo=math.floor(pos);hi=math.ceil(pos);return v[lo]+(v[hi]-v[lo])*(pos-lo)
def fmt(v):return f'{v:g}'
def summary(v):
    v=[x for x in v if x is not None]
    return f'{fmt(quantile(v,.5))} ({fmt(quantile(v,.25))}-{fmt(quantile(v,.75))}); n={len(v)}' if v else 'Unavailable'
paths=json.loads((ROOT/'inputs/source_paths.json').read_text())
md=read(ANALYSIS/'included_samples_descriptive113.csv')
fixed=[r for r in read(paths['table1_clinical']) if r['group'] != 'D']
assert len({r['Sample'] for r in fixed}) == len(fixed)
records={r['Sample']:r for r in fixed}
assert set(records) == {r['Sample'] for r in md if r['group'] != 'D'}
assert all(records[r['Sample']]['group'] == r['group'] for r in md if r['group'] != 'D')
from table1 import build_table1
results=build_table1(read(paths['table1_clinical']))
with (ROOT/'Table1.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['Characteristic','MV ARDS','Nonintubated ARDS','Pneumonia without ARDS','Healthy controls']);w.writerows(results)
# Supplementary table exports.
files={
'TableS1_DE_summary.csv':'gene_summary.csv','TableS2_Hallmark_all_contrasts.csv':'hallmark_all_contrasts.csv',
'TableS3a_module_sizes.csv':'wgcna/module_sizes.csv',
'TableS3b_all_patient_overlap.csv':'wgcna/all_patient_focal_overlap.csv','TableS3c_module_effects.csv':'wgcna/module_group_effects.csv',
'TableS4a_clinical_associations.csv':'clinical/clinical_associations.csv','TableS4b_platelet_adjustment.csv':'clinical/platelet_count_adjustment.csv',
'TableS5a_external_score.csv':'wgcna/GSE66890_focal_score_effect.csv','TableS5b_external_kME.csv':'wgcna/GSE66890_gene_module_correspondence.csv'}
for dst,src in files.items():shutil.copyfile(ANALYSIS/src,TABLES/dst)
# Keep submission tables consistent even when reading older nine-contrast caches.
for name in ['TableS1_DE_summary.csv', 'TableS2_Hallmark_all_contrasts.csv']:
    path = TABLES/name
    records = read(path)
    fields = list(records[0])
    records = [row for row in records if row['contrast'] != 'BC_D']
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
(ROOT/'Table1_notes.md').write_text('Continuous variables are median (IQR), with available n. Categorical variables are n (%). Clinical association analyses use the discovery cohort with complete cases for each variable.\n')
