import sys, json, os, re, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET, subprocess
# usage: python3 harvest.py queries.txt [max]   (one arXiv API query per line)
R='/Users/adriansager/playground/riemann'; D=R+'/.scratch/conjfals'
qf=sys.argv[1]; n=int(sys.argv[2]) if len(sys.argv)>2 else 40
store=D+'/harvest.json'
H=json.load(open(store)) if os.path.exists(store) else {}
out=subprocess.run(['rg','-o','--no-filename',r'\d{4}\.\d{4,5}','Riemann','NOGO.md','TREE.md','REFUTATIONS.md','DISCOVERIES.md','COROLLARIES.md','RECORDS.md','discoveries'],capture_output=True,text=True,cwd=R).stdout
seen=set(out.split())
ns={'a':'http://www.w3.org/2005/Atom'}
for q in [l.strip() for l in open(qf) if l.strip() and not l.startswith('#')]:
    url='http://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':q,'sortBy':'submittedDate','sortOrder':'descending','max_results':n})
    data=None
    for k in range(3):
        try:
            data=urllib.request.urlopen(url,timeout=60).read(); break
        except Exception as e:
            print('  retry',q[:50],str(e)[:60]); time.sleep(15)
    if data is None: print('FAIL',q); continue
    new=0; tot=0
    for e in ET.fromstring(data).findall('a:entry',ns):
        i=re.sub(r'v\d+$','',e.find('a:id',ns).text.rsplit('/',1)[-1])
        d=e.find('a:published',ns).text[:10]
        if d<'2023': continue
        tot+=1
        if i in H: H[i]['q'].append(q); continue
        t=' '.join(e.find('a:title',ns).text.split())
        a=' '.join(e.find('a:summary',ns).text.split())
        H[i]={'t':t,'d':d,'a':a,'q':[q],'seen':i in seen}
        new+=1
    print(f'Q {q[:90]} -> {tot} (new {new})')
    time.sleep(3.5)
json.dump(H,open(store,'w'))
print('total',len(H))
