# PYAudit: the two counterexamples to the S_3 chain of Patel-Yang Lemma 2.3 (TeX l. 490-499) and the 2609.00537 (S3-bound) remark.
# uv run --no-project --with mpmath python3 s3_examples.py
import mpmath as mp
mp.mp.dps = 30
def chain(al, be):
    nus = [n for n in range(int(mp.floor(al)) - 2, int(mp.ceil(be)) + 3) if al + 0.5 < n < be - 0.5]
    s_lo = sum(1/(n - al) for n in nus); s_hi = sum(1/(be - n) for n in nus)
    printed_n = sum(1/(n - mp.mpf(1)/2) for n in range(1, int(mp.ceil(be - al - 0.5)) + 1) if n < be - al - 0.5)
    psi = mp.digamma(be - al) - mp.digamma(mp.mpf(1)/2)
    S3 = (s_lo + s_hi)/mp.pi
    print('alpha=%s beta=%s nus=%s: sum 1/(nu-alpha)=%s | printed n-sum=%s | psi(b-a)-psi(1/2)=%s | |S3|=%s <= (2/pi)log(b-a)+2=%s' % (
        al, be, nus, mp.nstr(s_lo, 8), mp.nstr(printed_n, 8), mp.nstr(psi, 8), mp.nstr(S3, 8), mp.nstr(2/mp.pi*mp.log(be - al) + 2, 8)))
chain(mp.mpf('0.4'), mp.mpf('3.6'))
chain(mp.mpf('0.495'), mp.mpf('3.505'))
chain(mp.mpf('0.4999999'), mp.mpf('1.5000001'))
# 2609.00537 Cor. (Explicit_B_estimate) setting: 0 < f'(b) < 1, nu = 1 .. floor(f'(a)) - 1, S_3 as in PY (2.32)
al, be = mp.mpf('0.99'), mp.mpf('5.5')
S3 = sum(1/(be - n) + 1/(n - al) for n in range(1, int(mp.floor(be))))/mp.pi
print("2609.00537 (S3-bound): f'(b)=0.99, f'(a)=5.5: S_3 = %s vs (2/pi)log(f'(a)-f'(b)) + 1.251 = %s" % (mp.nstr(S3, 8), mp.nstr(2/mp.pi*mp.log(be - al) + mp.mpf('1.251'), 8)))
al, be = mp.mpf('0.4999999'), mp.mpf('1.5000001')
print("2609.00537-type constant at beta-alpha -> 1+ (PY range): |S3| -> 4/pi = %s > 1.251" % mp.nstr(4/mp.pi, 8))
# exact sup over beta-alpha = D in (1, inf) of [true sup |S_3| - printed (2/pi)log D - 2]: for D in (M, M+1] the true sup is
# (2/pi)(psi(M+1/2) - psi(1/2)) (both sums -> psi(M+1/2)-psi(1/2) as nu_1-alpha, beta-nu_M -> 1/2+), the printed bound increases in D,
# so the sup is approached at D -> M+.
vals = [(2/mp.pi*(mp.digamma(M + mp.mpf(1)/2) - mp.digamma(mp.mpf(1)/2)) - 2/mp.pi*mp.log(M) - 2, M) for M in range(1, 10001)]
print('sup_D [true - printed] = %s (approached at D -> %d+); limit M -> inf: %s' % (mp.nstr(max(vals)[0], 8), max(vals)[1], mp.nstr(2/mp.pi*(mp.euler + 2*mp.log(2)) - 2, 8)))
print('printed (2/pi)log D + 2 < 0 for D < e^{-pi} = %s, where S_3 is empty (= 0)' % mp.nstr(mp.e**(-mp.pi), 6))
