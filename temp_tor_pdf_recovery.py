#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, io, json, re, subprocess, sys, time, unicodedata
from pathlib import Path
import requests
from pypdf import PdfReader
import urllib3
urllib3.disable_warnings()

OUT=Path('notes_master_math_france_tor');OUT.mkdir(parents=True,exist_ok=True)
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126 Safari/537.36'
PROXIES={'http':'socks5h://127.0.0.1:9050','https':'socks5h://127.0.0.1:9050'}
RECORDS=[
{'num':1,'title':'Cours de probabilités','path':'01 - Proba-stat/01.1 Probabilites M1/01 - Cours de probabilites - Yves Coudene (2015).pdf','urls':['https://perso.lpsm.paris/~coudene/probabilites.pdf','https://www.lpsm.paris/pageperso/coudene/probabilites.pdf','http://www.proba.jussieu.fr/pageperso/coudene/probabilites.pdf'],'pages':[107],'keywords':['cours de probabil','yves coudene']},
{'num':19,'title':'Modélisation et statistique bayésienne computationnelle','path':'01 - Proba-stat/01.5 Bayesien et MCMC/19 - Modelisation et statistique bayesienne computationnelle - Nicolas Bousquet (2026).pdf','urls':['https://perso.lpsm.paris/~bousquet/poly-complet-2026-V1.pdf','https://www.lpsm.paris/pageperso/bousquet/poly-complet-2026-V1.pdf'],'pages':[142],'keywords':['bayes','bousquet']},
{'num':33,'title':'Contrôle optimal : théorie et applications','path':'02 - Analyse, optimisation et outils/Optimisation et controle/33 - Controle optimal theorie et applications - Emmanuel Trelat.pdf','urls':['https://www.ljll.fr/trelat/enseignement/controlSU/livreopt2.pdf','https://www.ljll.math.upmc.fr/trelat/fichiers/livreopt2.pdf','https://www.ljll.fr/~trelat/fichiers/livreopt.pdf'],'pages':[270,263,250,246],'keywords':['controle optimal','emmanuel trelat']},
{'num':34,'title':'Méthodes mathématiques et numériques pour les plasmas','path':'03 - EDP et calcul scientifique/34 - Methodes mathematiques et numeriques pour les plasmas - Bruno Despres (2021).pdf','urls':['https://www.ljll.fr/~despres/BD_fichiers/m2_plasma.pdf','https://www.ljll.fr/despres/BD_fichiers/m2_plasma.pdf','https://www.ljll.math.upmc.fr/despres/BD_fichiers/m2_plasma.pdf'],'pages':[61],'keywords':['plasma','despres']},
{'num':35,'title':'Équations aux dérivées partielles elliptiques','path':'03 - EDP et calcul scientifique/35 - Equations aux derivees partielles elliptiques - Herve Le Dret (2010).pdf','urls':['https://www.ljll.fr/ledret/M2Elliptique/chapitre4.pdf','https://www.ljll.math.upmc.fr/ledret/M2Elliptique/chapitre4.pdf'],'pages':[28],'keywords':['galerkin','le dret']},
]

def norm(s):
 s=unicodedata.normalize('NFKD',s);s=''.join(c for c in s if not unicodedata.combining(c));return re.sub(r'\s+',' ',s.lower())

def new_identity():
 subprocess.run(['sudo','killall','-HUP','tor'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 time.sleep(8)

def validate(data,rec):
 pos=data.find(b'%PDF-')
 if pos<0 or pos>4096:raise ValueError(f'non PDF {data[:40]!r}')
 data=data[pos:];reader=PdfReader(io.BytesIO(data),strict=False);pages=len(reader.pages)
 if pages not in rec['pages']:raise ValueError(f'pages inattendues {pages}')
 text=norm(' '.join((p.extract_text() or '') for p in reader.pages[:min(8,pages)]))
 if not any(k in text for k in rec['keywords']):raise ValueError(f'mots-cles absents {text[:300]!r}')
 return data,pages,text[:1000]

def download(url):
 r=requests.get(url,headers={'User-Agent':UA,'Accept':'application/pdf,*/*;q=0.8'},proxies=PROXIES,timeout=(30,240),verify=False,allow_redirects=True)
 r.raise_for_status();return r.content,r.url

def recover(rec):
 errors=[]
 for cycle in range(8):
  if cycle:new_identity()
  for url in rec['urls']:
   try:
    raw,final=download(url);data,pages,preview=validate(raw,rec);d=OUT/rec['path'];d.parent.mkdir(parents=True,exist_ok=True);d.write_bytes(data)
    return {**rec,'status':'ok','source_url':final,'source_kind':f'tor-cycle-{cycle}','actual_pages':pages,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'preview':preview,'detail':''}
   except Exception as e:errors.append(f'cycle {cycle} {url}: {type(e).__name__}: {e}')
 return {**rec,'status':'error','source_url':'','source_kind':'tor','actual_pages':0,'bytes':0,'sha256':'','preview':'','detail':' | '.join(errors)[-20000:]}

results=[]
for rec in RECORDS:
 x=recover(rec);results.append(x);print(f"[{x['num']:02d}] {x['status']} pages={x['actual_pages']} bytes={x['bytes']} {x['title']}",flush=True)
fields=['num','title','path','urls','pages','keywords','status','source_url','source_kind','actual_pages','bytes','sha256','preview','detail'];rows=[]
for x in results:
 y=dict(x);y['urls']=' | '.join(y['urls']);y['pages']=' | '.join(map(str,y['pages']));y['keywords']=' | '.join(y['keywords']);rows.append(y)
with (OUT/'manifest.csv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
(OUT/'manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
err=[x for x in results if x['status']!='ok'];print(f'SUCCES={len(results)-len(err)} ECHECS={len(err)}');sys.exit(1 if err else 0)
