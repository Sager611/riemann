import json, urllib.request, re, html, sys, time
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
for rid in sys.argv[1:]:
    for k in range(3):
        try:
            d=json.load(urllib.request.urlopen(urllib.request.Request('https://zenodo.org/api/records/'+rid,headers=UA),timeout=60)); break
        except Exception as e:
            print('retry',e); time.sleep(6)
    md=d['metadata']
    desc=' '.join(re.sub(r'<[^>]+>',' ',html.unescape(md.get('description',''))).split())
    print('==',rid,md.get('publication_date'),md.get('title'))
    print('creators:',[c.get('name') for c in md.get('creators',[])][:4])
    print(desc[:2000])
    print('files', [(f['key'], f['size']) for f in d.get('files',[])])
    time.sleep(2)
