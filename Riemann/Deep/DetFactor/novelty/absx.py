import sys, urllib.request, re, html, time
UA={'User-Agent':'Mozilla/5.0 (novelty scan; riemann repo)'}
for i in sys.argv[1:]:
    t=None
    for k in range(3):
        try:
            t=urllib.request.urlopen(urllib.request.Request('https://arxiv.org/abs/'+i,headers=UA),timeout=60).read().decode('utf-8','replace'); break
        except Exception as e:
            print('retry',i,e); time.sleep(8)
    if not t: continue
    ti=re.search(r'<h1 class="title mathjax"><span class="descriptor">Title:</span>(.*?)</h1>',t,re.S)
    ab=re.search(r'<blockquote class="abstract mathjax">\s*<span class="descriptor">Abstract:</span>(.*?)</blockquote>',t,re.S)
    au=re.search(r'<div class="authors">(.*?)</div>',t,re.S)
    da=re.search(r'\[Submitted on ([^\]]*)\]',t)
    clean=lambda s: ' '.join(html.unescape(re.sub(r'<[^>]+>','',s)).split()) if s else ''
    print('==',i,'|',clean(ti.group(1)) if ti else '', '|', clean(au.group(1))[:120] if au else '', '|', da.group(1) if da else '')
    print('  ',clean(ab.group(1))[:1500] if ab else 'NO ABSTRACT')
    time.sleep(2)
