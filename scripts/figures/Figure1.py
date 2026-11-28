from pathlib import Path
import re, csv, collections
from html import escape
p=Path('results/figures/Figure1')
p.mkdir(parents=True,exist_ok=True)
defs=Path('scripts/assets/figure1_defs.svg').read_text()
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="114.4mm" viewBox="0 100 1200 750" role="img" aria-labelledby="title desc"><title id="title">Plasma cfRNA profiling across ARDS clinical states</title><desc id="desc">Four color-coded clinical groups feed into blood collection, plasma isolation, cfRNA extraction, library preparation and sequencing, followed by pathway and co-expression analyses.</desc><defs>{defs}</defs><rect width="1200" height="920" fill="white"/>']
def text(x,y,t,size=18,bold=False,anchor='middle',color='#243849'):
 s.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" style="font-family:Arial,Helvetica,sans-serif;font-size:{size}px;fill:{color};font-weight:{700 if bold else 400}">{escape(t)}</text>')
def box(x,y,w,h,fill='#FFFFFF',stroke='#D6E0E5',r=10):s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="1.3"/>')
def use(n,x,y,z=44,c='#547E91'):s.append(f'<use href="#{n}" x="{x}" y="{y}" width="{z}" height="{z}" style="color:{c}"/>')
def path(d,arrow=False):s.append(f'<path d="{d}" class="arrow"'+(' marker-end="url(#arrow)"' if arrow else '')+'/>')
# Cohort row and laboratory row have independent, balanced spacing.
with Path('results/inputs/analysis_samples.csv').open() as f:
 counts=collections.Counter(r['group'] for r in csv.DictReader(f))
assert set(counts) == set('ABCD'), 'Figure 1 requires prepared samples in groups A/B/C/D'
patient_n = sum(counts[g] for g in 'ABC')
total_n = sum(counts.values())
rows=[(165,'#B38A56',['Mechanically ventilated','ARDS'],str(counts['A'])),(455,'#93749E',['Non-mechanically','ventilated ARDS'],str(counts['B'])),(745,'#527D9A',['Pneumonia','without ARDS'],str(counts['C'])),(1035,'#8B999F',['Healthy controls'],str(counts['D']))]
for x,c,lines,n in rows:
 use('person',x-20,107,40,c)
 for i,l in enumerate(lines):text(x,173+i*22,l,18)
 text(x,226,'n = '+n,18,True,color=c)
path('M165 246 V259 H1035 V246')
path('M600 259 V286',True)
steps=[('blood',115,['Blood','collection']),('plasma',309,['Plasma','isolation']),('rna',503,['cfRNA','extraction']),('library',697,['Library','preparation']),('sequencer',891,['RNA','sequencing'])]
for name,x,labels in steps:
 use(name,x-26,298,52)
 text(x,378,labels[0],17);text(x,400,labels[1],17)
for x in [115,309,503,697,891]:path(f'M{x+43} 324 H{x+151}',True)
for i in range(5):
 for j in range(5):
  c=['#E4EDF1','#B8CDD7','#7FA5B7','#517F95'][(i+j*2)%4]
  s.append(f'<rect x="{1060+j*11}" y="{298+i*11}" width="9" height="9" fill="{c}"/>')
text(1085,378,'Expression',17);text(1085,400,'profiles',17)
# Broad-to-focused analysis, with healthy controls used descriptively.
path('M1085 410 V423 H600 V439',True)
text(600,453,f'Disease background: patients and healthy controls (n = {total_n})',21,True)
text(600,483,f'Patient-only analysis (n = {patient_n}): overall ARDS vs pneumonia without ARDS',19)
path('M600 496 V518',True)
text(600,544,'Comparisons across the three patient groups',18,True)
path('M600 558 V576 M315 576 H885 M315 576 V593',True);path('M885 576 V593',True)
box(35,601,560,113,'#F5F9FC','#C5D7E2');box(615,601,550,113,'#F9F6FB','#D9CCE0')
use('pathway',56,617,34,'#527D9A');text(104,642,'Pathway analysis',22,True,'start')
text(315,682,'Hallmark and expanded biological programs',18)
use('network',635,617,34,'#93749E');text(685,642,'Co-expression discovery',22,True,'start')
text(890,682,'Modules · Stability · Functional annotation',18)
path('M315 714 V735 H885 V714');path('M600 735 V749')
path('M215 749 H985')
for x in [215,600,985]:path(f'M{x} 749 V763',True)
for x,t in [(35,'Clinical-state patterns'),(420,'Clinical associations'),(805,'Independent dataset')]:
 box(x,771,360,57,'#FAFBFC','#D8E1E6')
 # Matching line icons: clinical states, associations, external replication.
 s.append(f'<g transform="translate({x+19},785)" fill="none" stroke="#547E91" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">')
 if x == 35:
  s.append('<circle cx="7" cy="5" r="3"/><circle cx="24" cy="5" r="3"/><path d="M2 23V16Q2 11 7 11Q12 11 12 16V23M19 23V16Q19 11 24 11Q29 11 29 16V23M13 7H18M15 5L18 7L15 9"/>')
 elif x == 420:
  s.append('<path d="M2 1V25H31M7 20L25 5"/><circle cx="8" cy="17" r="1.7"/><circle cx="14" cy="15" r="1.7"/><circle cx="20" cy="8" r="1.7"/><circle cx="27" cy="7" r="1.7"/>')
 else:
  s.append('<path d="M3 2H20V26H3Z M7 7H16M7 12H13M7 17H12"/><circle cx="24" cy="19" r="8" fill="#FAFBFC"/><path d="M20 19L23 22L28 16"/>')
 s.append('</g>')
 text(x+62,806,t,21,True,'start')
s.append('</svg>')
(p/'Figure1.svg').write_text('\n'.join(s))

# Export every Figure 1 format in this script; R/magick is shared with the other figures.
import subprocess
subprocess.run(['Rscript','-e',r"""
suppressPackageStartupMessages(library(magick))
f1 <- 'results/figures/Figure1/Figure1.svg'
im <- image_read(f1,density=300)
image_write(im,sub('[.]svg$','.png',f1),format='png')
image_write(image_read(f1,density=600),sub('[.]svg$','.tiff',f1),format='tiff',compression='lzw')
image_write(im,sub('[.]svg$','.pdf',f1),format='pdf')
"""],check=True)
