
import json,re,os
ROOT='/Users/adriansager/playground/riemann'
d=json.load(open(ROOT+'/.scratch/oaisext/closure.json'))
down=set(d['down']); closure=set(d['closure'])
def p(m): return os.path.join(ROOT,*m.split('.'))+'.lean'
names=['NaturalState','plainEnergy','energy','positiveSlotRow','rowTwistedSum','idealRowHom','FixedInducingRow','Profiles','diagonalControl','control','ZeroAt','CertifiedBand','certified_terminal','terminal_certificate','ZeroBound','Radial','RadialData','internalQ','actualSextic','zeroMaximum','familyZeros']
defre=re.compile(r'^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(def|structure|abbrev|class|inductive|theorem|lemma|instance)\s+([^\s:(\[{]+)',re.M)
where={}
for m in closure:
    src=open(p(m),encoding='utf-8').read()
    for mm in defre.finditer(src):
        nm=mm.group(2)
        base=nm.split('.')[-1]
        if base in names:
            where.setdefault(base,[]).append((mm.group(1),m.replace('Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.',''),'DOWN' if m in down else 'up'))
for k in names:
    print(k, where.get(k))
# defs in down files mentioning beta
print('--- defs/structures in down files whose text block mentions beta ---')
wb=re.compile(r'\bbeta\b')
for m in sorted(down):
    src=open(p(m),encoding='utf-8').read()
    # split into top-level declaration blocks
    blocks=re.split(r'\n(?=(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(?:def|structure|abbrev|class|inductive|theorem|lemma|instance)\s)',src)
    for b in blocks:
        mm=defre.match(b.lstrip('\n'))
        if not mm: continue
        if mm.group(1) in ('def','structure','abbrev','class','inductive','instance') and wb.search(b):
            print(m.replace('Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.',''), mm.group(1), mm.group(2))

