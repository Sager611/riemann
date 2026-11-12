# gl2audit: float kappa(R) = max over lam, Q of 1 - log M/R at R outside and across the certified box (sanity; myceil Ritz)
import sys
args = sys.argv[1:]; sys.argv = ['m', args[0], args[1], args[2]]
import myceil as MY
mp = MY.mp
for Rs in args[3].split(','):
    R = mp.mpf(Rs); F, lam, bb = MY.Malt(R, lam0=mp.mpf('0.8'), it=200)
    print(f"R={Rs}: kappa={mp.nstr(1-mp.log(F)/R, 11)} lam*={mp.nstr(lam, 6)}", flush=True)
