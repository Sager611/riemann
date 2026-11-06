import re, sys, json, urllib.request, urllib.parse, time
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
def get(url):
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=60) as r:
                return r.read().decode('utf-8','replace')
        except Exception as e:
            err=str(e); time.sleep(3)
    return 'FETCHFAIL '+err
mode=sys.argv[1]
if mode=='api':
    t=get('https://hexagonmath.org/api/v1'); print(t[:3000])
elif mode=='search':
    for q in sys.argv[2:]:
        t=get('https://hexagonmath.org/search?q='+urllib.parse.quote(q))
        ids=sorted(set(re.findall(r'href="/(\d{4}\.\d{5}v\d+)"',t)))
        # titles: try to capture link text
        pairs=re.findall(r'href="/(\d{4}\.\d{5}v\d+)"[^>]*>(.*?)</a>',t,re.S)
        seen=set(); out=[]
        for i,ti in pairs:
            ti=re.sub(r'<[^>]+>','',ti); ti=' '.join(ti.split())
            if i in seen or not ti: continue
            seen.add(i); out.append(i+' '+ti[:110])
        print('Q',q,'->',len(ids),'ids'); print('\n'.join(out[:25]))
