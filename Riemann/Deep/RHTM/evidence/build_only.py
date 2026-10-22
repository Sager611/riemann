import sys, json, hashlib
from pathlib import Path
import importlib
# Run from the directory holding the extracted appendix modules (see NOTE.md).
ROOT=Path.cwd()
sys.path.insert(0,str(ROOT/'vendor')); sys.path.insert(0,str(ROOT))
def _m(name): return importlib.import_module(name)
RHBuilder=_m('rh120_source').RHBuilder
save,count=_m('builder').save,_m('builder').count
minimize=_m('minimize').minimize
canonical,read,write=_m('tables').canonical,_m('tables').read,_m('tables').write
_w=_m('window_certificate'); generate,verify,verify_projection=_w.generate,_w.verify,_w.verify_projection
optimize,quot=_m('compatible_merge').optimize,_m('compatible_merge').quot

def compile_named(b):
    raw=b.build(False);sm,mp,rounds=minimize(raw)
    for (q,s),(w,d,n) in raw.items():
        assert sm[mp[q],s]==(w,d,mp.get(n,'halt'))
    tab,ids,missing=canonical(sm);unmerged,rawids,_=canonical(raw)
    mapping=[ids[mp[q]] for q,i in sorted(rawids.items(),key=lambda x:x[1])]
    for (q,s),(w,d,n) in unmerged.items():
        assert tab[mapping[q],s]==(w,d,-1 if n<0 else mapping[n])
    return raw,tab,{q:ids[mp[q]] for q in mp},unmerged,mapping,rounds,missing

b=RHBuilder();raw,tab,names,unmerged,bmap,rounds,missing=compile_named(b)
print('raw states',count(raw),'canonical predecessor states',len(tab)//2,'missing',missing)
cert=generate(tab,2);info=verify(tab,cert)
groups,trace=optimize(tab,info['reachable_read_overapproximation'])
final,projection=quot(tab,info['reachable_read_overapproximation'],groups)
verify_projection(tab,final,cert,projection)
print('final',len(final), 'groups', str(groups)[:300])
paper=read(ROOT/'machines/RH_120_paper.tm')
print('final == paper table:', final==paper)
save(raw,ROOT/'machines/RH_121_raw_named.tm')
write(tab,ROOT/'machines/RH_121_predecessor.tm')
write(unmerged,ROOT/'machines/RH_122_unmerged.tm')
write(final,ROOT/'machines/RH_120.tm')
json.dump(names,open(ROOT/'machines/RH_121.names.json','w'),indent=1)
json.dump({q:projection[i] for q,i in names.items()},open(ROOT/'machines/RH_120.names.json','w'),indent=1)
json.dump(b.dispatch_certificate,open(ROOT/'machines/RH_121.dispatch.json','w'),indent=1)
print('window info', {k:v for k,v in info.items() if k!='reachable_read_overapproximation'})
print('root', b.root)
