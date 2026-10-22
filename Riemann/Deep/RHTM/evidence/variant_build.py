import sys, json
from pathlib import Path
ROOT=Path.cwd()  # directory holding the extracted appendix modules (see NOTE.md)
sys.path.insert(0,str(ROOT/'vendor')); sys.path.insert(0,str(ROOT))
import importlib
_src=importlib.import_module('rh120_source')
RHBuilder, subroutine = _src.RHBuilder, _src.subroutine
minimize=importlib.import_module('minimize').minimize
_tables=importlib.import_module('tables')
canonical, write = _tables.canonical, _tables.write
class Inflate(RHBuilder):
    @subroutine
    def seed_deficit(self):
        # FALSIFIED VARIANT: f = 2(n-1)+1 instead of n-1 (deficit inflated by n)
        return [self.f.inc, self.while_decnz(self.j,[self.v.inc,self.f.inc,self.f.inc])]
b=Inflate(); raw=b.build(False); sm,mp,rounds=minimize(raw); tab,ids,missing=canonical(sm)
write(tab,Path('variant_inflate.tm'))
names={q:ids[mp[q]] for q in mp}
roles={'3.m.inc':'leaf00','3.m.decnz':'leaf01','3.j.inc':'leaf10','3.j.decnz':'leaf11','3.x.inc':'leaf20','3.x.decnz':'leaf21','3.f.inc':'leaf30','3.f.decnz':'leaf31','3.v.inc':'leaf40','3.v.decnz':'leaf41','2c.reg.return_0_1':'ret0','2c.reg.return_1_1':'ret1','5.root.2':'root2','N%d'%b.root:'rootN'}
with open('roles_inflate.txt','w') as fh:
    for k,r in roles.items(): fh.write('%s %d\n'%(r,names[k]))
print('variant states',len(tab)//2,'root',b.root,'missing',missing)
