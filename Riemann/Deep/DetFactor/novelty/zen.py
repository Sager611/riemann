import sys, json, urllib.request, urllib.parse, re, html, time
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
def get(u):
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=60) as r: return r.read().decode('utf-8','replace')
        except Exception as e: err=str(e); time.sleep(5)
    return None
mode=sys.argv[1]
if mode=='q':
    for q in sys.argv[2:]:
        t=get('https://zenodo.org/api/records?'+urllib.parse.urlencode({'q':q,'sort':'mostrecent','size':25}))
        if t is None: print('FAIL',q); continue
        d=json.loads(t)
        print('Q',q,'total',d['hits']['total'])
        for h in d['hits']['hits']:
            md=h['metadata']; print('  ',h['id'],md.get('publication_date'),md.get('title','')[:110])
elif mode=='rec':
    pat=re.compile(sys.argv[2],re.I)
    for rid in sys.argv[3:]:
        t=get('https://zenodo.org/api/records/'+rid)
        if t is None: print('FAIL',rid); continue
        d=json.loads(t); md=d['metadata']
        desc=re.sub(r'<[^>]+>',' ',html.unescape(md.get('description','')))
        desc=' '.join(desc.split())
        print('==',rid,md.get('publication_date'),md.get('title','')[:100],'| files:',[f['key'] for f in d.get('files',[])][:5])
        for m in pat.finditer(desc):
            s=max(0,m.start()-150); print('   ...',desc[s:m.end()+150])
