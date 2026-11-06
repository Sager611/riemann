import re, sys, json, urllib.request, urllib.parse, time, os, io, zipfile, tarfile, gzip
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
def getb(url):
    err=''
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=90) as r:
                return r.read(), r.headers.get('Content-Type')
        except Exception as e:
            err=str(e); time.sleep(3)
    return None, err
for ident in sys.argv[1:]:
    b,ct=getb('https://hexagonmath.org/api/v1/works/%s/source'%ident)
    os.makedirs('hxsrc/'+ident,exist_ok=True)
    if b is None: print(ident,'FAIL',ct); continue
    print(ident, ct, len(b))
    if ct and 'json' in ct:
        d=json.loads(b)
        # tree listing: try to find files
        s=json.dumps(d)[:800]; print(s)
        open('hxsrc/%s/tree.json'%ident,'w').write(json.dumps(d))
    else:
        open('hxsrc/%s/blob'%ident,'wb').write(b)
