import random
random.seed(20261119)
def sqf(u):
    p=2
    while p*p<=u:
        if u%(p*p)==0: return False
        p+=1
    return True
def fund(d):
    r=d%16
    if r in (3,7,11,15): return sqf(d)
    if r in (4,8): return sqf(d//4)
    return False
out=[]
while len(out)<300:
    d=random.randrange(10**5,10**9)
    if fund(d): out.append(d)
open('small_d.txt','w').write('\n'.join(map(str,out))+'\n')
