import sys, os, re, json, time, gzip, tarfile, io, subprocess
# usage: python3 dl_scan.py ids.txt   -> downloads e-prints (3.2 s apart) into tex/<id>/, scans, writes scan.json
D='/Users/adriansager/playground/riemann/.scratch/conjfals'
T=D+'/tex'; os.makedirs(T,exist_ok=True)
ids=[l.split()[0] for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
def fetch(i):
    d=f'{T}/{i}'
    if os.path.isdir(d) and os.listdir(d): return 'cached'
    os.makedirs(d,exist_ok=True)
    f=f'{T}/{i}.src'
    for k in range(2):
        r=subprocess.run(['curl','-s','-L','-A','Mozilla/5.0 (riemann conjfals)','-o',f,'-w','%{http_code}',f'https://arxiv.org/e-print/{i}'],capture_output=True,text=True,timeout=120)
        if r.stdout.strip()=='200': break
        time.sleep(15)
    if not os.path.exists(f): return 'nofile'
    b=open(f,'rb').read(); os.remove(f)
    if b[:4]==b'%PDF': return 'pdf-only'
    try: b=gzip.decompress(b)
    except Exception: pass
    try:
        tf=tarfile.open(fileobj=io.BytesIO(b))
        n=0
        for m in tf.getmembers():
            if m.isfile() and m.name.lower().endswith(('.tex','.bbl','.ltx')) and m.size<5e6:
                open(f'{d}/'+re.sub(r'[^\w.-]','_',m.name),'wb').write(tf.extractfile(m).read()); n+=1
        return f'tar {n}'
    except Exception:
        if b[:4]==b'%PDF': return 'pdf-only'
        open(f'{d}/main.tex','wb').write(b); return 'single'
S=json.load(open(D+'/scan.json')) if os.path.exists(D+'/scan.json') else {}
PH=re.compile(r'we conjecture|is conjectured|conjecture that|numerical (evidence|computations?|experiments?) (suggest|indicate|support)|we expect that|open (problem|question)|it seems (that|likely)|we believe|we ask whether|it would be interesting to (know|determine|decide)|remains open',re.I)
for k,i in enumerate(ids):
    st=fetch(i)
    if st not in ('cached',): time.sleep(3.2)
    d=f'{T}/{i}'
    envnames={'conjecture','conj','question','problem','openproblem','oproblem','prob','ques','qn','conjec','con','hypothesis','observation'}
    texs=[x for x in sorted(os.listdir(d)) if x.endswith(('.tex','.ltx'))]
    txt={}
    for x in texs:
        try: txt[x]=open(f'{d}/{x}',encoding='utf-8',errors='replace').read()
        except Exception: pass
    for x,s in txt.items():
        for m in re.finditer(r'\\newtheorem\*?\{(\w+)\}(?:\[\w+\])?\{([^}]*)\}',s):
            if re.search(r'conj|question|problem|hypothes|observ|claim|speculat',m.group(2),re.I): envnames.add(m.group(1).lower())
        for m in re.finditer(r'\\(?:declaretheorem|newmdtheoremenv)(?:\[[^\]]*name=([^,\]]*)[^\]]*\])?\{(\w+)\}',s):
            if m.group(1) and re.search(r'conj|question|problem',m.group(1),re.I) or re.search(r'conj|quest|prob',m.group(2),re.I): envnames.add(m.group(2).lower())
    envs=[]; phr=[]
    for x,s in txt.items():
        for m in re.finditer(r'\\begin\{(\w+\*?)\}(.*?)\\end\{\1\}',s,re.S):
            if m.group(1).rstrip('*').lower() in envnames:
                ln=s.count('\n',0,m.start())+1
                envs.append([x,ln,m.group(1),' '.join(m.group(2).split())[:900]])
        for ln,line in enumerate(s.split('\n'),1):
            if PH.search(line) and not line.lstrip().startswith('%'):
                phr.append([x,ln,' '.join(line.split())[:300]])
    S[i]={'st':st,'envs':envs,'phr':phr[:40],'nenv':len(envs),'nphr':len(phr)}
    print(f'{k+1}/{len(ids)} {i} {st} envs={len(envs)} phr={len(phr)}',flush=True)
    json.dump(S,open(D+'/scan.json','w'))
