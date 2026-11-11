import os, re, json, sys
D='/Users/adriansager/playground/riemann/.scratch/conjfals'; T=D+'/tex'
PH=re.compile(r'we conjecture|is conjectured|conjecture that|numerical (evidence|computations?|experiments?) (suggest|indicate|support)|we expect that|open (problem|question)|it seems (that|likely)|we believe|we ask whether|it would be interesting to (know|determine|decide)|remains open',re.I)
S={}
for i in sorted(os.listdir(T)):
    d=f'{T}/{i}'
    if not os.path.isdir(d): continue
    envnames={'conjecture','conj','question','problem','openproblem','oproblem','prob','ques','qn','conjec','con','observation','conjecturex','prob*'}
    txt={}
    for x in sorted(os.listdir(d)):
        if x.endswith(('.tex','.ltx')):
            txt[x]=open(f'{d}/{x}',encoding='utf-8',errors='replace').read()
    for x,s in txt.items():
        for m in re.finditer(r'\\newtheorem\*?\{(\w+)\}(?:\[\w+\])?\{([^}]*)\}',s):
            if re.search(r'conj|question|problem|observ|speculat',m.group(2),re.I): envnames.add(m.group(1).lower())
        for m in re.finditer(r'\\declaretheorem(?:\[([^\]]*)\])?\{(\w+)\}',s):
            if re.search(r'conj|quest|prob',(m.group(1) or '')+m.group(2),re.I): envnames.add(m.group(2).lower())
    envs=[]; phr=[]
    for x,s in txt.items():
        for nm in envnames:
            for m in re.finditer(r'\\begin\{'+re.escape(nm)+r'\*?\}(.*?)\\end\{'+re.escape(nm)+r'\*?\}',s,re.S|re.I):
                ln=s.count('\n',0,m.start())+1
                envs.append([x,ln,nm,' '.join(m.group(1).split())[:1200]])
        for ln,line in enumerate(s.split('\n'),1):
            if PH.search(line) and not line.lstrip().startswith('%'):
                phr.append([x,ln,' '.join(line.split())[:300]])
    envs.sort(key=lambda e:(e[0],e[1]))
    S[i]={'envs':envs,'phr':phr[:60],'nenv':len(envs),'nphr':len(phr)}
json.dump(S,open(D+'/scan2.json','w'))
for i,v in S.items():
    if v['nenv'] or v['nphr']: print(i,'envs',v['nenv'],'phr',v['nphr'])
