import json,re,os,difflib
ROOT='/Users/adriansager/playground/riemann'
man=json.load(open(ROOT+'/.scratch/oaisext/fork_manifest.json'))
def p(m): return os.path.join(ROOT,*m.split('.'))+'.lean'
bad=[]; changed_lines=0; files_changed=0; out=[]
for e in man:
    o=open(p(e['orig']),encoding='utf-8').read(); f=open(p(e['fork']),encoding='utf-8').read()
    assert 'betaEff' not in o
    fl=f.split('\n'); assert fl[2].startswith('-- oaisext fork'); f2='\n'.join(fl[:2]+fl[3:])
    f2=re.sub(r'^(\s*import\s+)Riemann\.Deep\.OAISext\.',r'\1Riemann.Deep.QRH.',f2,flags=re.M)
    f2=f2.replace('\nimport Riemann.Deep.QRH.BetaEff','')
    f2=f2.replace('HeckeZeroSupremum.betaEff','HeckeZeroSupremum.beta')
    f2=re.sub(r'(?<![\w.])betaEff\b','beta',f2)
    if f2!=o:
        d=list(difflib.unified_diff(o.split('\n'),f2.split('\n'),lineterm='',n=0))
        bad.append((e['orig'],d))
    # raw line diff counts orig vs fork
    dd=[l for l in difflib.unified_diff(o.split('\n'),f.split('\n'),lineterm='',n=0) if (l.startswith('+') or l.startswith('-')) and not l.startswith('+++') and not l.startswith('---')]
    if dd: files_changed+=1; changed_lines+=sum(1 for l in dd if l.startswith('+'))
print('files',len(man),'files differing from original (incl. header/imports):',files_changed,'added/changed lines:',changed_lines)
print('files whose diff is NOT pure (header + import rewrite + beta->betaEff):',len(bad))
for m,d in bad:
    print('==',m)
    for l in d: print('   ',l)
