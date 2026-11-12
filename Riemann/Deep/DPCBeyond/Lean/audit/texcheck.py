
import tarfile, gzip, io, os, re, sys
d='.scratch/dpcblean/tex'
os.makedirs(d, exist_ok=True)
def load(path):
    data=open(path,'rb').read()
    try:
        t=tarfile.open(fileobj=io.BytesIO(data))
        texs={m.name: t.extractfile(m).read().decode('utf8','replace') for m in t.getmembers() if m.name.endswith('.tex')}
        return texs
    except Exception:
        try:
            return {'main.tex': gzip.decompress(data).decode('utf8','replace')}
        except Exception:
            return {'main.tex': data.decode('utf8','replace')}
def show(name, texs, ranges, kw):
    # pick the biggest tex
    fn=max(texs, key=lambda k: len(texs[k]))
    L=texs[fn].split('\n')
    print(f'## {name}: {fn}, {len(L)} lines')
    for a,b in ranges:
        hits=[]
        for i in range(a-1, min(b, len(L))):
            s=L[i].strip()
            if any(k.lower() in s.lower() for k in kw) :
                hits.append(f'{i+1}: {s[:110]}')
        print(f'[{a}-{b}] ' + (' | '.join(h for h in hits[:3]) if hits else '(no keyword) first: ' + L[a-1].strip()[:100]))
lz=load('.scratch/dpcblean/lz.tgz'); dpc=load('.scratch/dpcblean/dpc.tgz')
kw=['begin{lemma','begin{proposition','begin{theorem','label','N(\\sigma','Fujii','zero','density','Riesz','Rodgers','explicit','Lemma','Prop','Hölder','Holder','mollif','\\theta','sign','E_T','covering','consecutive','phi(1']
show('Lamzouri 2611.14486', lz, [(176,177),(221,230),(240,250),(251,284)], kw)
show('DPC 2610.28527v2', dpc, [(262,270),(335,345),(614,626),(641,643),(673,683),(686,702),(743,756),(753,753),(760,772),(774,787)], kw)

