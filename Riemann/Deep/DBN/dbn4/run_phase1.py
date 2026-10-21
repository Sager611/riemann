
import json, sys, time, gzip, hashlib, shutil
from pathlib import Path
from fractions import Fraction as F
CERT=Path('/Users/adriansager/playground/riemann/.scratch/dbn3/src/anc/certificate')
sys.path.insert(0,str(CERT/'code'))
import produce_classical_portfolio as prod  # pyright: ignore[reportMissingImports]
import audit_classical_portfolio as aud  # pyright: ignore[reportMissingImports]
here=Path('.').resolve()
cat=json.load(open('catalogue.json'))
ev=here/'evidence'; ev.mkdir(exist_ok=True)
names=['SOURCE_OLD_ENHANCED.json','SOURCE_RETUNED.json']
new=[{'r':k['r'],'rho':k['rho'],'lambda':k['lam'],'self_cap':k['self_cap'],'dbn4_certificate':True,'leaves':k['leaves']} for k in cat['kernels']]
for n in names:
    src=json.loads((CERT/'inputs/refined'/n).read_text())
    orig=json.loads((CERT/'inputs/refined'/n).read_text())
    assert len(src['kernels'])==57
    src['kernels']=src['kernels']+new
    # everything except the kernel list is byte-for-byte the refined Planat source
    for key in orig:
        if key!='kernels': assert src[key]==orig[key],key
    assert src['kernels'][:57]==orig['kernels']
    (ev/n).write_text(json.dumps(src,separators=(',',':'),sort_keys=True)+'\n')
print('kernels per provider',57+len(new))
t=time.time()
out=here/'histories'; out.mkdir(exist_ok=True)
res=prod.produce(ev,out,[0,1],57+len(new),32,'classical_dbn4')
print('produced',res['upper_decimal'],res['steps'],'%.1fs'%(time.time()-t),flush=True)
paths=[ev/n for n in names]; S=[json.loads(p.read_text()) for p in paths]; hashes=[aud.sha(p) for p in paths]
with gzip.open(out/'classical_dbn4.json.gz','rt') as f: H=json.load(f)
t=time.time(); A=aud.audit(H,S,hashes); print('AUDIT',A['PASS'],A['upper_exact'],A['upper_decimal'],A['steps'],'%.1fs'%(time.time()-t))
json.dump(A,open(out/'classical_dbn4_AUDIT.json','w'),indent=2,sort_keys=True)
neg=aud.negative_controls(H,S,hashes); print('NEGATIVE_CONTROLS',neg['PASS'],len(neg['cases']))
used=sorted({r['kernel'] for r in H['records']})
print('distinct kernels used',len(used),'new used',sum(1 for k in used if k>=57))
