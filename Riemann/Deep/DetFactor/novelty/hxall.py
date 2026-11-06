import re, sys, json, urllib.request, urllib.parse, time, os
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
def get(url):
    err=''
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=60) as r:
                return r.read().decode('utf-8','replace')
        except Exception as e:
            err=str(e); time.sleep(4)
    return None
out='hx_works.json'
works=[]
cur=None
for it in range(80):
    url='https://hexagonmath.org/api/v1/works?limit=100' + ('&cursor='+urllib.parse.quote(cur) if cur else '')
    t=get(url)
    if t is None: print('FAIL at',it); break
    d=json.loads(t)
    for w in d['works']:
        works.append({k:w.get(k) for k in ('identifier','title','abstract','subject','publishedAt','tags')})
    pg=d['page']
    if not pg.get('has_more'): break
    cur=pg.get('next_cursor')
json.dump(works,open(out,'w'))
print('works',len(works), 'total', d['page'].get('total'))
pat=re.compile(sys.argv[1],re.I)
for w in works:
    s=(w['title'] or '')+' '+(w['abstract'] or '')
    if pat.search(s):
        print(w['identifier'], w['publishedAt'][:10] if w['publishedAt'] else '', '|', (w['title'] or '')[:120])
