import urllib.request, re, html, time, sys
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
pat=re.compile(sys.argv[1],re.I)
for mon in sys.argv[2:]:
    for skip in (0,250,500):
        url='https://arxiv.org/list/cs.SC/%s?skip=%d&show=250'%(mon,skip)
        t=None
        for k in range(3):
            try:
                t=urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=60).read().decode('utf-8','replace'); break
            except Exception as e:
                print('retry',e); time.sleep(10)
        if not t: print('FAIL',mon,skip); continue
        ids=re.findall(r'<dt>.*?href\s*=\s*"/abs/(\d{4}\.\d{4,5})',t,re.S)
        titles=[html.unescape(re.sub(r'\s+',' ',x)).strip() for x in re.findall(r"list-title[^>]*>\s*<span class=.descriptor.>Title:</span>(.*?)</div>",t,re.S)]
        tot=re.search(r'Total of (\d+) entries',t)
        print(mon,skip,'n',len(ids),'total',tot.group(1) if tot else '?')
        for i,ti in zip(ids,titles):
            if pat.search(ti): print('  ',i,ti[:140])
        if not tot or int(tot.group(1))<=skip+250: break
        time.sleep(3)
