#!/usr/bin/env python3
from __future__ import annotations
import io, json, re, sys
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup
from pypdf import PdfReader
import urllib3
urllib3.disable_warnings()

OUT=Path('trelat_dspace');OUT.mkdir(exist_ok=True)
UA='Mozilla/5.0 Chrome/126 Safari/537.36'
H={'User-Agent':UA,'Accept':'text/html,application/pdf,*/*;q=0.8'}
BASES=['https://di.univ-blida.dz','http://di.univ-blida.dz']
HANDLE_PATHS=['/jspui/handle/123456789/38391','/xmlui/handle/123456789/38391','/handle/123456789/38391']
CANDIDATES=[]
for base in BASES:
 for p in HANDLE_PATHS:CANDIDATES.append(urljoin(base,p))
 for p in [
 '/jspui/bitstream/123456789/38391/1/2-515-317.pdf',
 '/jspui/bitstream/123456789/38391/2/2-515-317.pdf',
 '/jspui/bitstream/123456789/38391/3/2-515-317.pdf',
 '/xmlui/bitstream/handle/123456789/38391/2-515-317.pdf?sequence=1&isAllowed=y',
 '/xmlui/bitstream/handle/123456789/38391/2-515-317.pdf?sequence=2&isAllowed=y',
 '/bitstream/123456789/38391/1/2-515-317.pdf',
 '/server/api/core/items/38391/content',
 ]:CANDIDATES.append(urljoin(base,p))

def validate(data):
 pos=data.find(b'%PDF-')
 if pos<0 or pos>4096:raise ValueError(f'non PDF {data[:40]!r}')
 data=data[pos:];r=PdfReader(io.BytesIO(data),strict=False);pages=len(r.pages)
 text=' '.join((p.extract_text() or '') for p in r.pages[:5]).lower()
 if pages<200 or 'contr' not in text or 'optimal' not in text:raise ValueError(f'validation title/pages failed pages={pages} text={text[:200]!r}')
 return data,pages,text[:1000]

s=requests.Session();s.headers.update(H);errors=[];seen=set();queue=list(CANDIDATES)
while queue:
 url=queue.pop(0)
 if url in seen:continue
 seen.add(url)
 try:
  r=s.get(url,timeout=(15,120),verify=False,allow_redirects=True)
  errors.append(f'{url} -> {r.status_code} {r.url} {r.headers.get("content-type")} {len(r.content)}')
  if r.status_code>=400:continue
  try:
   data,pages,preview=validate(r.content)
   path=OUT/'33 - Controle optimal theorie et applications - Emmanuel Trelat.pdf';path.write_bytes(data)
   (OUT/'result.json').write_text(json.dumps({'status':'ok','url':r.url,'pages':pages,'bytes':len(data),'preview':preview},ensure_ascii=False,indent=2),encoding='utf-8')
   print('SUCCESS',r.url,pages,len(data));sys.exit(0)
  except Exception:pass
  if 'html' in (r.headers.get('content-type') or '').lower() or r.content.lstrip().startswith(b'<'):
   soup=BeautifulSoup(r.content,'html.parser')
   for a in soup.find_all('a',href=True):
    href=a.get('href','')
    text=a.get_text(' ',strip=True)
    if '2-515-317' in href or '2-515-317' in text or 'bitstream' in href.lower() or 'download' in href.lower():
     queue.append(urljoin(r.url,href))
 except Exception as e:errors.append(f'ERROR {url}: {type(e).__name__}: {e}')
(OUT/'result.json').write_text(json.dumps({'status':'error','tried':len(seen),'errors':errors},ensure_ascii=False,indent=2),encoding='utf-8')
print('\n'.join(errors));sys.exit(1)
