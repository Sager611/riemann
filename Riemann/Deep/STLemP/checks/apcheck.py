
def ap(p):
    s=0
    for x in range(p):
        v=(4*x*x*x-4*x+1)%p
        if v==0: continue
        s+= 1 if pow(v,(p-1)//2,p)==1 else -1
    return -s
print([ (p,ap(p)) for p in [5,7,11,13,17,19,23,29,31,41,43]])
