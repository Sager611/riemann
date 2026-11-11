import json,re,os,sys
ROOT='/Users/adriansager/playground/riemann'
QRH='Riemann.Deep.QRH.'; FORK='Riemann.Deep.OAISext.'
d=json.load(open(ROOT+'/.scratch/oaisext/closure.json'))
down=set(d['down'])
def p(m): return os.path.join(ROOT,*m.split('.'))+'.lean'
HDR='-- oaisext fork (Riemann/Deep/OAISext/REPORT.md): HeckeZeroSupremum.beta replaced by HeckeZeroSupremum.betaEff (max of beta and 51/100); imports of forked modules point to Riemann.Deep.OAISext.'
manifest=[]; warn=[]
tot_q=0; tot_bare=0
for m in sorted(down):
    src=open(p(m),encoding='utf-8').read()
    lines=src.split('\n')
    assert lines[0].startswith('-- Ported') and lines[1].startswith('-- Apache'), m
    s=src
    def rimp(mm):
        mod=mm.group(2)
        if mod in down: return mm.group(1)+FORK+mod[len(QRH):]
        return mm.group(0)
    s=re.sub(r'^(\s*import\s+)(\S+)',rimp,s,flags=re.M)
    nfork_imp=sum(1 for i in re.findall(r'^\s*import\s+(\S+)',s,flags=re.M) if i.startswith(FORK))
    pat_q=re.compile(r'(?<![\w.])HeckeZeroSupremum\.beta\b')
    hyp_ok=len(re.findall(r'\(51/100:ℝ\)≤HeckeZeroSupremum\.beta\b',s))+len(re.findall(r'2\*HeckeZeroSupremum\.beta-1≤κ',s))
    s,nq=pat_q.subn('HeckeZeroSupremum.betaEff',s)
    nb=0
    opens=any(re.match(r'\s*open\b.*\bHeckeZeroSupremum\b',l) for l in lines)
    if opens:
        s,nb=re.subn(r'(?<![\w.])beta\b','betaEff',s)
    tot_q+=nq; tot_bare+=nb
    if m.endswith('Moments.PrimeGlobal'):
        a='import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Hecke.ZeroSupremum'
        assert a in s
        s=s.replace(a,a+'\nimport Riemann.Deep.OAISext.BetaEff')
        old='    exact zero_re_le_beta (χ i) (by linarith [hz.1]) (Or.inr (hχ i)) hz.2.2.2'
        assert old in s
        s=s.replace(old,'    exact (zero_re_le_beta (χ i) (by linarith [hz.1]) (Or.inr (hχ i)) hz.2.2.2).trans beta_le_betaEff')
        assert 'betaEff i hbeta beta_le_one' in s
        s=s.replace('betaEff i hbeta beta_le_one','betaEff i hbeta betaEff_le_one')
    ls_=s.split('\n'); s='\n'.join(ls_[:2]+[HDR]+ls_[2:])
    dst=os.path.join(ROOT,'Riemann','Deep','OAISext',*m[len(QRH):].split('.'))+'.lean'
    os.makedirs(os.path.dirname(dst),exist_ok=True)
    open(dst,'w',encoding='utf-8').write(s)
    manifest.append({'orig':m,'fork':FORK+m[len(QRH):],'qualified':nq,'hyp_pattern':hyp_ok,'bare':nb,'fork_imports':nfork_imp})
    if hyp_ok!=nq: warn.append((m,nq,hyp_ok))
json.dump(manifest,open(ROOT+'/.scratch/oaisext/fork_manifest.json','w'),indent=0)
print('forked',len(manifest),'qualified',tot_q,'bare',tot_bare,'files with substitutions',sum(1 for x in manifest if x['qualified'] or x['bare']),'warn',warn)
