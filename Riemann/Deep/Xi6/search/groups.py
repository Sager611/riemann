#!/usr/bin/env python3
# groups.py: write multiplication tables for small groups
import sys, itertools
def write(fn, n, mul):
    with open(fn,'w') as f:
        f.write(f"{n}\n")
        for a in range(n): f.write(' '.join(map(str,mul[a]))+'\n')
def cyclic(n): return n, [[(a+b)%n for b in range(n)] for a in range(n)]
def abelian(ns):
    els=list(itertools.product(*[range(x) for x in ns])); idx={e:i for i,e in enumerate(els)}
    mul=[[idx[tuple((x+y)%q for x,y,q in zip(a,b,ns))] for b in els] for a in els]
    return len(els), mul
def dihedral(k):  # order 2k: (r,s) r in Z_k, s in {0,1}: element = rot^r * ref^s
    els=[(r,s) for s in range(2) for r in range(k)]; idx={e:i for i,e in enumerate(els)}
    def mm(a,b):
        r1,s1=a; r2,s2=b
        return ((r1 + (r2 if s1==0 else -r2))%k, (s1+s2)%2)
    return len(els), [[idx[mm(a,b)] for b in els] for a in els]
def semidirect(p, a, d):  # Z_p x| Z_d, (x,i)*(y,j) = (x + a^i y, i+j), a of order d mod p
    els=[(x,i) for i in range(d) for x in range(p)]; idx={e:t for t,e in enumerate(els)}
    def mm(u,v):
        x,i=u; y,j=v
        return ((x + pow(a,i,p)*y)%p, (i+j)%d)
    return len(els), [[idx[mm(u,v)] for v in els] for u in els]
def gf_add(p, e):  # additive group of GF(p^e) = Z_p^e
    return abelian([p]*e)
if __name__=='__main__':
    kind=sys.argv[1]; fn=sys.argv[-1]
    if kind=='cyc': n,mul=cyclic(int(sys.argv[2]))
    elif kind=='ab': n,mul=abelian([int(x) for x in sys.argv[2].split('x')])
    elif kind=='dih': n,mul=dihedral(int(sys.argv[2]))
    elif kind=='sd': n,mul=semidirect(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
    else: raise SystemExit('kind?')
    assert mul[0]==list(range(n))
    write(fn,n,mul)

