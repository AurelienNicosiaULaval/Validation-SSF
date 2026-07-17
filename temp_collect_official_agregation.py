#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, io, json, re, shutil, time
from collections import deque
from pathlib import Path
from urllib.parse import urljoin, urlparse, urlunparse

import requests
from bs4 import BeautifulSoup
from pypdf import PdfReader

OUT=Path('agregation_externe_officiel_complet')
OUT.mkdir(parents=True,exist_ok=True)
STARTS=[
 ('Jury actuel - archives','https://agreg.org/index.php?id=archives',3),
 ('Jury actuel - modelisation','https://agreg.org/index.php?id=modelisation',3),
 ('Ancien jury - sujets','https://old.agreg.org/sujets.html',3),
 ('Ancien jury - accueil','https://old.agreg.org/',3),
 ('Ministere - programme 2027','https://www.devenirenseignant.gouv.fr/les-programmes-des-concours-d-enseignants-du-second-degre-de-la-session-2027-1662',2),
 ('Ministere - sujets 2026','https://www.devenirenseignant.gouv.fr/sujets-et-rapports-des-jurys-concours-de-l-agregation-de-2026-1659',2),
]
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126 Safari/537.36'
S=requests.Session();S.headers.update({'User-Agent':UA,'Accept':'text/html,application/pdf,*/*;q=0.8','Accept-Language':'fr-FR,fr;q=0.9'})
last={}
records=[];hashes={};seen=set()
PDF_HINT=re.compile(r'\.pdf(?:$|[?#])|/data/uploads/|sujet|rapport|model|mod\d|mg\d|ap\d|programme',re.I)
HTML_HINT=re.compile(r'archive|rapport|sujet|model|session|agreg|concours',re.I)

def norm(u):
 p=urlparse(u);return urlunparse((p.scheme.lower(),p.netloc.lower(),p.path,p.params,p.query,''))

def polite_get(url,stream=False):
 host=urlparse(url).netloc.lower();wait=1.6-(time.monotonic()-last.get(host,0))
 if wait>0:time.sleep(wait)
 for attempt in range(8):
  try:
   r=S.get(url,timeout=(15,120),allow_redirects=True,stream=stream)
   last[host]=time.monotonic()
   if r.status_code==429:
    retry=r.headers.get('Retry-After')
    delay=min(60,float(retry)) if retry and retry.isdigit() else min(60,4*(attempt+1))
    r.close();time.sleep(delay);continue
   if r.status_code>=500:
    r.close();time.sleep(min(30,3*(attempt+1)));continue
   return r
  except Exception:
   if attempt==7:raise
   time.sleep(min(30,3*(attempt+1)))
 raise RuntimeError('échec après reprises')

def clean(s):
 s=re.sub(r'[\\/:*?"<>|\x00-\x1f]+',' - ',s or '')
 s=re.sub(r'\s+',' ',s).strip(' .-_')
 return (s or 'document')[:180]

def validate_pdf(data):
 pos=data.find(b'%PDF-')
 if pos<0 or pos>4096:raise ValueError('signature PDF absente')
 data=data[pos:];reader=PdfReader(io.BytesIO(data),strict=False);n=len(reader.pages)
 if n<1:raise ValueError('aucune page')
 title=''
 try:title=str((reader.metadata or {}).get('/Title') or '').strip()
 except Exception:pass
 return data,n,title

def save_pdf(label,source_page,url):
 try:
  r=polite_get(url,stream=True)
  if r.status_code>=400:
   records.append({'source':label,'source_page':source_page,'url':url,'status':f'http_{r.status_code}','filename':'','bytes':0,'pages':0,'sha256':'','detail':r.reason or ''});r.close();return
  content=bytearray()
  for chunk in r.iter_content(262144):
   if chunk:content.extend(chunk)
   if len(content)>150*1024*1024:raise ValueError('fichier >150 Mio')
  final=r.url;cd=r.headers.get('content-disposition','');ctype=r.headers.get('content-type','');r.close()
  raw=bytes(content)
  raw,pages,title=validate_pdf(raw)
  digest=hashlib.sha256(raw).hexdigest()
  if digest in hashes:
   records.append({'source':label,'source_page':source_page,'url':url,'status':'duplicate','filename':'','bytes':len(raw),'pages':pages,'sha256':digest,'detail':hashes[digest]});return
  name=Path(urlparse(final).path).name or Path(urlparse(url).path).name or 'document.pdf'
  m=re.search(r'filename\*=UTF-8\'\'([^;]+)',cd,re.I)
  if not m:m=re.search(r'filename="?([^";]+)',cd,re.I)
  if m:name=m.group(1)
  if not name.lower().endswith('.pdf'):name=(title or name)+'.pdf'
  name=clean(name[:-4])+'.pdf'
  folder=OUT/clean(label);folder.mkdir(parents=True,exist_ok=True)
  dest=folder/name
  if dest.exists():dest=folder/f'{dest.stem} [{digest[:8]}].pdf'
  dest.write_bytes(raw);rel=str(dest.relative_to(OUT));hashes[digest]=rel
  records.append({'source':label,'source_page':source_page,'url':url,'status':'ok','filename':rel,'bytes':len(raw),'pages':pages,'sha256':digest,'detail':ctype})
  print('PDF',pages,rel,flush=True)
 except Exception as e:
  records.append({'source':label,'source_page':source_page,'url':url,'status':'error','filename':'','bytes':0,'pages':0,'sha256':'','detail':f'{type(e).__name__}: {e}'})

def crawl(label,start,max_depth):
 q=deque([(norm(start),0)]);host=urlparse(start).netloc.lower();pages=0
 while q and pages<1000:
  url,depth=q.popleft()
  if url in seen or depth>max_depth:continue
  seen.add(url)
  try:
   r=polite_get(url)
   if r.status_code>=400:
    records.append({'source':label,'source_page':'','url':url,'status':f'http_{r.status_code}','filename':'','bytes':0,'pages':0,'sha256':'','detail':r.reason or ''});r.close();continue
   ctype=r.headers.get('content-type','').lower();data=r.content;final=norm(r.url);r.close()
   if data.find(b'%PDF-') in range(0,4097):save_pdf(label,url,url);continue
   if 'html' not in ctype and not data.lstrip().startswith((b'<html',b'<!DOCTYPE')):continue
   pages+=1;soup=BeautifulSoup(data.decode('utf-8','replace'),'lxml')
   for a in soup.find_all('a',href=True):
    href=a.get('href','').strip();text=' '.join(a.get_text(' ',strip=True).split())
    if not href or href.startswith(('#','mailto:','javascript:')):continue
    target=norm(urljoin(final,href));pt=urlparse(target)
    if PDF_HINT.search(target) or target.lower().endswith('.pdf'):
     if target not in seen:
      seen.add(target);save_pdf(label,final,target)
    elif pt.netloc.lower()==host and depth<max_depth and HTML_HINT.search(target+' '+text):
     q.append((target,depth+1))
  print('PAGE',pages,label,final,flush=True)
  except Exception as e:
   records.append({'source':label,'source_page':'','url':url,'status':'error','filename':'','bytes':0,'pages':0,'sha256':'','detail':f'{type(e).__name__}: {e}'})

for label,url,d in STARTS:crawl(label,url,d)
fields=['source','source_page','url','status','filename','bytes','pages','sha256','detail']
with (OUT/'MANIFEST.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
(OUT/'MANIFEST.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
summary={'records':len(records),'pdfs':sum(r['status']=='ok' for r in records),'duplicates':sum(r['status']=='duplicate' for r in records),'errors':sum(r['status'].startswith('http_') or r['status']=='error' for r in records),'pages':sum(int(r['pages'] or 0) for r in records if r['status']=='ok'),'bytes':sum(int(r['bytes'] or 0) for r in records if r['status']=='ok')}
(OUT/'SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'README.md').write_text('# Archives officielles de l’agrégation externe de mathématiques\n\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
shutil.make_archive(OUT.name,'zip',OUT)
print(json.dumps(summary,ensure_ascii=False),flush=True)
