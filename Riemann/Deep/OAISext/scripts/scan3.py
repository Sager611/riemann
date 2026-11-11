
import json,re,os
ROOT='/Users/adriansager/playground/riemann'
d=json.load(open(ROOT+'/.scratch/oaisext/closure.json'))
down=sorted(d['down'])
def p(m): return os.path.join(ROOT,*m.split('.'))+'.lean'
for m in down:
    src=open(p(m),encoding='utf-8').read()
    lines=src.split('\n')
    hits=[]
    for i,l in enumerate(lines):
        for o in re.finditer(r'(?<![\w.])beta\b',l):
            hits.append((i+1,l.strip()[:150]))
    if hits:
        opens=[l.strip() for l in lines if re.match(r'\s*open\b',l) and 'HeckeZeroSupremum' in l]
        print('==',m.replace('Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.',''),'opens HZS:',bool(opens))
        for h in hits: print('  ',h)

