import urllib.request, urllib.parse, json, sys, time
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
for q in sys.argv[1:]:
    u='https://api.archives-ouvertes.fr/search/?'+urllib.parse.urlencode({'q':q,'fq':'producedDateY_i:[2025 TO 2026]','fl':'halId_s,title_s,producedDate_s','rows':20,'sort':'producedDate_s desc'})
    try:
        d=json.load(urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=60))
        r=d['response']; print('Q',q,'n',r['numFound'])
        for x in r['docs']: print('  ',x.get('halId_s'),x.get('producedDate_s'),(x.get('title_s') or [''])[0][:120])
    except Exception as e: print('FAIL',q,e)
    time.sleep(1)
