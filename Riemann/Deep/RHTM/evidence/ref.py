"""Independent op-level reference interpreter of the paper'self five-register source
(rh120_source = rh_source + width regrouping), written from the source semantics,
not from the builder. Emits the primitive-op trace hash used by tmsim.c."""
import sys
M64=(1<<64)-1
m,j,x,f,v=0,1,2,3,4
class Halt(Exception): pass
class RM:
    def __init__(self, regs=(0,0,0,0,0), variant=None):
        self.r=list(regs); self.h=0; self.ops=0; self.variant=variant; self.marks={}
    def inc(self,a):
        self.h=(self.h*1000003 + a*4) & M64; self.ops+=1; self.r[a]+=1
    def dec(self,a):
        if self.r[a]>0:
            self.h=(self.h*1000003 + a*4+2) & M64; self.ops+=1; self.r[a]-=1; return True
        self.h=(self.h*1000003 + a*4+3) & M64; self.ops+=1; return False
    def mark(self,name): self.marks.setdefault(name,(self.ops,tuple(self.r)))
    # macros
    def move(self,a,b):
        while self.dec(a): self.inc(b)
    def fold(self):
        while self.dec(j):
            self.inc(v); self.inc(x)
        self.move(v,j)
    def transform(self):
        self.move(x,m)
        while True:
            self.fold()
            if self.dec(m): self.inc(x); continue
            break
    def divide(self):
        self.move(j,v); self.move(x,m)
        while self.dec(m):
            if self.dec(v): self.inc(j)
            else: self.move(j,v); self.inc(x)
    def product_step(self):
        self.inc(j); self.divide()
        if self.dec(v): self.inc(v); self.transform()
        self.transform(); self.transform()
    def seed(self):
        if self.variant=='inflate':
            self.inc(f)
            while self.dec(j): self.inc(v); self.inc(f); self.inc(f)
        else:
            while self.dec(j): self.inc(v); self.inc(f)
    def digits(self):
        while self.dec(x):
            self.divide(); self.dec(f)
        self.move(v,j)
    def halve(self,a):
        while self.dec(a):
            if a==j: self.inc(m)
            if self.dec(a): self.inc(v)
        self.move(v,a)
    def width(self):
        while self.dec(j):
            self.inc(m); self.halve(j); self.halve(j); self.halve(f); self.inc(x)
        self.mark('compare')
        while self.dec(x): self.dec(f)
        if self.dec(f): raise Halt()
    def stage(self):
        self.mark('stage')
        self.move(m,f)
        while True:
            self.product_step()
            if not self.dec(f): break
        self.mark('seed'); self.seed()
        self.mark('digits'); self.digits()
        self.mark('width'); self.width()

def run_stages(maxn, variant=None):
    rm=RM(variant=variant); out=[]
    gops=0
    for n in range(2,maxn+1):
        rm.h=0; rm.ops=0; rm.marks={}
        try:
            rm.stage(); halted=False
        except Halt:
            halted=True
        out.append(dict(n=n,ops=rm.ops,hash='%016x'%rm.h,regs=tuple(rm.r),halted=halted,
                        marks={k:(gops+o,r) for k,(o,r) in rm.marks.items()}))
        gops+=rm.ops
        if halted: break
    return out

if __name__=='__main__':
    maxn=int(sys.argv[1]); variant=sys.argv[2] if len(sys.argv)>2 else None
    for d in run_stages(maxn,variant):
        print(d['n'],d['ops'],d['hash'],d['regs'],'HALT' if d['halted'] else '',d['marks'])
