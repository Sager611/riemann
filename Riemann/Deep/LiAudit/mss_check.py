# Compare lambda_tilde = lambda_n - lambda_bar_n (FFT tables) with the values in Merlini-Sala-Sala 1903.06233 Sec. 2 table.
import sys
def rd(f):
    d = {}
    for l in open(f):
        if l[0] != '#': p = l.split(); d[int(p[0])] = float(p[1])
    return d
lam, bar = rd(sys.argv[1]), rd(sys.argv[2])
mss = {5: 1.458268, 10: 1.323802, 20: 0.670652, 31: 0.733544, 100: 0.628752, 500: 2.663502, 1000: 1.756264, 2000: 10.76850,
       3000: -2.09000, 200: 3.600, 760: -6.33, 1870: -11.226, 3300: -9.900, 5080: 16.830, 5500: 16.500, 6500: 13.000, 8000: 16.000}
for n in sorted(mss): print('n=%5d ours %+.6f  MSS %+.6f' % (n, lam[n] - bar[n], mss[n]))
