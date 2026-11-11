import mpmath as mp, sys
mp.mp.dps=20
a,b=int(sys.argv[1]),int(sys.argv[2])
with open('zeros_%d_%d.txt'%(a,b),'w') as f:
    for n in range(a,b+1):
        f.write('%d %s\n'%(n, mp.nstr(mp.zetazero(n).imag,18))); f.flush()
