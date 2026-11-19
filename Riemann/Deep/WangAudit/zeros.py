# compute first N zeta zeros (ordinates) with mpmath; write .scratch/wangaudit/zeros.txt
import sys, time
import mpmath as mp
mp.mp.dps = 20
N = int(sys.argv[1]); out = sys.argv[2]
t0 = time.time()
with open(out, 'w') as fh:
    for n in range(1, N+1):
        z = mp.zetazero(n)
        fh.write(f"{n} {mp.nstr(z.imag, 18)}\n")
        if n % 250 == 0:
            fh.flush(); print(n, round(time.time()-t0,1), flush=True)
print('done', N, round(time.time()-t0,1))
