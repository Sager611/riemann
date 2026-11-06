import re, sys, json, urllib.request, urllib.parse, time
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
def get(url):
    err=''
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=60) as r:
                return r.read().decode('utf-8','replace')
        except Exception as e:
            err=str(e); time.sleep(3)
    return 'FETCHFAIL '+err
t=get('https://hexagonmath.org/api/v1/works')
try:
    d=json.loads(t)
except Exception:
    print(t[:500]); sys.exit()
print(type(d), list(d.keys())[:20] if isinstance(d,dict) else len(d))
if isinstance(d,dict):
    for k,v in d.items():
        if isinstance(v,list): print(k,len(v)); print(json.dumps(v[0])[:1500] if v else '')
        else: print(k, str(v)[:200])
