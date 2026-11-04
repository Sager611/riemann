
# Independent exact check of the projection-body counterexample (product of two 10-simplices in R^20)
# Method: Pi(P) = sum_F [-a_F/2, a_F/2] (a_F = area vector of facet F); vol(zonotope sum_i [-v_i/2,v_i/2]) = sum_{|S|=d} |det(v_S)|.
# Facet area vectors computed from scratch (no product formula used): for a simplex conv(p_0..p_d), the area vector of the facet opposite p_k
# is  -(1/(d-1)!) * cofactor-normal; we compute it as  a_k = - d * vol(simplex) * grad of barycentric coordinate lambda_k.
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, factorial
import sys

def det(M):
    # exact Gaussian elimination over Fractions
    M=[row[:] for row in M]; n=len(M); d=Fr(1)
    for c in range(n):
        p=next((r for r in range(c,n) if M[r][c]!=0),None)
        if p is None: return Fr(0)
        if p!=c: M[c],M[p]=M[p],M[c]; d=-d
        d*=M[c][c]; inv=1/M[c][c]
        for r in range(c+1,n):
            f=M[r][c]*inv
            if f: M[r]=[a-f*b for a,b in zip(M[r],M[c])]
    return d

def simplex_area_vectors(pts):
    # pts: d+1 points in R^d. Barycentric lambda_k(x) is affine; grad lambda_k = row k of inverse of [[p_j - p_0]] etc.
    d=len(pts)-1
    # matrix with columns (p_j - p_0), j=1..d
    A=[[Fr(pts[j][i]-pts[0][i]) for j in range(1,d+1)] for i in range(d)]
    # invert A
    n=d; aug=[A[i][:]+[Fr(int(i==j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p=next(r for r in range(c,n) if aug[r][c]!=0); aug[c],aug[p]=aug[p],aug[c]
        inv=1/aug[c][c]; aug[c]=[x*inv for x in aug[c]]
        for r in range(n):
            if r!=c and aug[r][c]!=0:
                f=aug[r][c]; aug[r]=[x-f*y for x,y in zip(aug[r],aug[c])]
    Ainv=[row[n:] for row in aug]
    vol=abs(det(A))/factorial(d)
    grads=[None]*(d+1)
    for k in range(1,d+1): grads[k]=Ainv[k-1]           # lambda_k = (A^{-1}(x-p0))_k
    grads[0]=[-sum(Ainv[k][i] for k in range(d)) for i in range(d)]
    # outward area vector of facet opposite p_k:  a_k = -d*vol*grad lambda_k
    return vol,[[ -d*vol*g for g in grads[k]] for k in range(d+1)]

def zonotope_volume(gens,d):
    return sum(abs(det([list(g) for g in S])) for S in combinations(gens,d))

def simplex(d): return [[0]*d]+[[int(i==j) for j in range(d)] for i in range(d)]

# sanity: check area vectors sum to zero (closed polytope) and simplex ratio formula for small d by brute force
for d in [2,3,4,5]:
    vol,av=simplex_area_vectors(simplex(d))
    assert all(sum(a[i] for a in av)==0 for i in range(d))
    R=zonotope_volume(av,d)/vol**(d-1)
    assert R==Fr((d+1)*d**d,factorial(d)), (d,R)
print('simplex ratio R_d=(d+1)d^d/d! verified by brute-force zonotope sums for d=2..5')
# also check product identity by brute force on small product T2 x T3 (R^5) using facets of the product built explicitly
def product_area_vectors(d1,d2):
    v1,a1=simplex_area_vectors(simplex(d1)); v2,a2=simplex_area_vectors(simplex(d2))
    return v1*v2,[a+[Fr(0)]*d2 for a in [[x*v2 for x in aa] for aa in a1]]+[[Fr(0)]*d1+[x*v1 for x in aa] for aa in a2]
for d1,d2 in [(2,2),(2,3),(3,3)]:
    vol,gens=product_area_vectors(d1,d2)
    R=zonotope_volume(gens,d1+d2)/vol**(d1+d2-1)
    R1=Fr((d1+1)*d1**d1,factorial(d1)); R2=Fr((d2+1)*d2**d2,factorial(d2))
    assert R==R1*R2,(d1,d2,R,R1*R2)
print('product identity R(AxB)=R(A)R(B) verified by brute force for (2,2),(2,3),(3,3)')
# main case d=20: 22 generators in R^20 -> 231 determinants of 20x20 rational matrices
vol,gens=product_area_vectors(10,10)
Z=zonotope_volume(gens,20)
R=Z/vol**19
c20=Fr(21*20**20,factorial(20))
q=R/c20
print('R_20(T10xT10)/c_20 =',q,'=',float(q))
print('matches 121*C(20,10)/(21*2^20):',q==Fr(121*comb(20,10),21*2**20),' and 22355476/22020096:',q==Fr(22355476,22020096),' >1:',q>1)
