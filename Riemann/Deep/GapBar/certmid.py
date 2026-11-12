"""gapbar: Arb certificate of h1 > 1 for a REAL polynomial P with the MIDPOINT mollifier M(s + i kappa pi/(2L)).
BM's Lemma main (l.370-384) holds uniformly for shifts << 1/L, including the mollifier shifts alpha3, beta3; with
alpha3 = -beta3 = i kappa pi/(2L) its kernel is BM's kernel times y^{-alpha3(x1+x2)-beta3(x3+x4)} = e^{ic(r2-r1)/2}
(r1 = X1+X2, r2 = X3+X4), the real twisted kernel of gtw.py.  So this is GapCeil's certr.py with gtw's Gamma arrays.
F(t) = e^{ivtL} zeta(s) zeta(s + i kappa pi/L) M(s + i kappa pi/(2L)) vanishes at t_n and t_n - kappa pi/L as BM need.
usage: certmid.py logs/optmid_....json kappa [D=80]"""
import certr
import gtw
certr.gamma_arrays_r = gtw.gamma_arrays_tw
if __name__ == '__main__':
    certr.main()
    print('  midpoint mollifier (twisted kernel); Arb Im parts of the twisted arrays contain 0:', gtw.IMCHECK)
