#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures as cf, csv, hashlib, io, json, subprocess, sys, tempfile
from pathlib import Path
from pypdf import PdfReader

OUT=Path('notes_master_math_france_lpsm_pageperso');OUT.mkdir(parents=True,exist_ok=True)
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126 Safari/537.36'
RECORDS=[
(1,'Cours de probabilités','01 - Proba-stat/01.1 Probabilites M1/01 - Cours de probabilites - Yves Coudene (2015).pdf',[
'https://www.lpsm.paris/pageperso/coudene/probabilites.pdf','http://www.lpsm.paris/pageperso/coudene/probabilites.pdf']),
(5,'Probabilités approfondies : martingales et chaînes de Markov','01 - Proba-stat/01.2 Processus, Markov et martingales/05 - Probabilites approfondies martingales et chaines de Markov - Thomas Duquesne (2012).pdf',[
'https://www.lpsm.paris/pageperso/broutin/teaching/4M011_poly_duquesne.pdf','https://www.lpsm.paris/pageperso/broutinn/teaching/4M011_poly_duquesne.pdf','http://www.lpsm.paris/pageperso/broutin/teaching/4M011_poly_duquesne.pdf']),
(11,'Statistique, Partie 2 : approche bayésienne','01 - Proba-stat/01.5 Bayesien et MCMC/11 - Statistique, Partie 2 approche bayesienne - Anna Ben-Hamou; Arnaud Guyader.pdf',[
'https://www.lpsm.paris/pageperso/aguyader/files/teaching/M1/PolycopiePartie2.pdf','http://www.lpsm.paris/pageperso/aguyader/files/teaching/M1/PolycopiePartie2.pdf']),
(15,'Calcul stochastique et processus de diffusion','01 - Proba-stat/01.3 Calcul stochastique et diffusions/15 - Calcul stochastique et processus de diffusion - Nicolas Fournier.pdf',[
'https://www.lpsm.paris/pageperso/nfournier/PolyCS.pdf','https://www.lpsm.paris/pageperso/fournier/PolyCS.pdf','http://www.lpsm.paris/pageperso/nfournier/PolyCS.pdf']),
(19,'Modélisation et statistique bayésienne computationnelle','01 - Proba-stat/01.5 Bayesien et MCMC/19 - Modelisation et statistique bayesienne computationnelle - Nicolas Bousquet (2026).pdf',[
'https://www.lpsm.paris/pageperso/bousquet/poly-complet-2026-V1.pdf','http://www.lpsm.paris/pageperso/bousquet/poly-complet-2026-V1.pdf']),
]

def validate(data):
 p=data.find(b'%PDF-')
 if p<0 or p>4096: raise ValueError(f'non PDF {data[:40]!r}')
 data=data[p:]; n=len(PdfReader(io.BytesIO(data),strict=False).pages)
 if n<1: raise ValueError('aucune page')
 return data,n

def curl(url):
 with tempfile.NamedTemporaryFile(delete=False) as f: tmp=Path(f.name)
 cmd=['curl','-L','--fail','--retry','3','--retry-all-errors','--connect-timeout','12','--max-time','100','-A',UA,'-H','Accept: application/pdf,*/*;q=0.8','-o',str(tmp),url]
 try:
  p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  if p.returncode: raise RuntimeError((p.stderr or p.stdout)[-1500:])
  return validate(tmp.read_bytes())
 finally: tmp.unlink(missing_ok=True)

def recover(rec):
 num,title,path,urls=rec;errors=[]
 for url in urls:
  try:
   data,pages=curl(url);d=OUT/path;d.parent.mkdir(parents=True,exist_ok=True);d.write_bytes(data)
   return {'num':num,'title':title,'path':path,'urls':urls,'status':'ok','source_url':url,'pages':pages,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'detail':''}
  except Exception as e:errors.append(f'{url}: {type(e).__name__}: {e}')
 return {'num':num,'title':title,'path':path,'urls':urls,'status':'error','source_url':'','pages':0,'bytes':0,'sha256':'','detail':' | '.join(errors)[-10000:]}

with cf.ThreadPoolExecutor(max_workers=5) as ex: results=list(ex.map(recover,RECORDS))
results.sort(key=lambda x:x['num'])
for x in results: print(f"[{x['num']:02d}] {x['status']} pages={x['pages']} bytes={x['bytes']} {x['title']}",flush=True)
fields=['num','title','path','urls','status','source_url','pages','bytes','sha256','detail'];rows=[]
for x in results:
 y=dict(x);y['urls']=' | '.join(y['urls']);rows.append(y)
with (OUT/'manifest.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
(OUT/'manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
err=[x for x in results if x['status']!='ok'];print(f'SUCCES={len(results)-len(err)} ECHECS={len(err)}');sys.exit(1 if err else 0)
