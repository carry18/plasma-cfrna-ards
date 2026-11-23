"""Parse the official GEO series matrix/platform files; no processed result imports."""
import csv, gzip, json
from pathlib import Path
out=Path('results/inputs')
paths=json.loads((out/'source_paths.json').read_text())
meta={};rows=[];inside=False
with gzip.open(paths['geo_matrix'],'rt') as f:
    for row in csv.reader(f,delimiter='\t'):
        if not row: continue
        if row[0].startswith('!series_matrix_table_begin'): inside=True;continue
        if row[0].startswith('!series_matrix_table_end'): break
        if inside: rows.append(row)
        elif row[0] in ['!Sample_title','!Sample_geo_accession']: meta[row[0]]=row[1:]
        elif row[0]=='!Sample_characteristics_ch1':
            key=row[1].split(': ',1)[0];meta[key]=[v.split(': ',1)[-1] for v in row[1:]]
with gzip.open(out/'GSE66890_expression.tsv.gz','wt') as f: csv.writer(f,delimiter='\t').writerows(rows)
with (out/'GSE66890_metadata.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(list(meta));w.writerows(zip(*meta.values()))
with gzip.open(paths['geo_platform'],'rt') as f:
    for line in f:
        if line.startswith('ID\t'): break
    else: raise ValueError('Missing platform header')
    mapping=[]
    for r in csv.DictReader(f,delimiter='\t',fieldnames=line.rstrip('\n').split('\t')):
        symbol=(r.get('Gene symbol') or '').strip()
        if symbol and symbol!='---' and not any(x in symbol for x in ['///','//',';','|',',']):mapping.append([r['ID'],symbol])
with (out/'GPL6244_unambiguous_mapping.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['probe','gene']);w.writerows(mapping)
