# Check the transported Jacobi matrix used in certify6.py (CJ): closed forms for diag(k) - V and off(k)^2, and that its
# eigenvalues are the roots of q_{U,V}(Y) = 2F1(-d, U; V; Y/U).  Not part of the proof.
import mpmath as mp
mp.mp.dps = 60
worst = 0
for (d, V, U) in [(6, 300, 600), (9, 500, 2000), (12, 1000, 1013), (7, 400, 10 ** 7), (10, 800, 1600)]:
    al, be = V - 1, U - V - d; S = al + be
    bk = [mp.mpf(be ** 2 - al ** 2) / ((2 * k + S) * (2 * k + S + 2)) for k in range(d)]
    ak2 = [mp.mpf(4 * k * (k + al) * (k + be) * (k + S)) / ((2 * k + S) ** 2 * ((2 * k + S) ** 2 - 1)) for k in range(1, d)]
    diag = [U * (1 - b) / 2 for b in bk]
    off = [U * mp.sqrt(x) / 2 for x in ak2]
    # closed forms
    for k in range(d):
        cf = mp.mpf(V * S * (d - 1) + 2 * k * (S + k + 1) * (U - 2 * V)) / ((2 * k + S) * (2 * k + S + 2))
        worst = max(worst, abs(diag[k] - V - cf))
    for k in range(1, d):
        cf2 = mp.mpf(U ** 2 * k * (k + V - 1) * (k + U - V - d) * (k + S)) / ((2 * k + S) ** 2 * ((2 * k + S) ** 2 - 1))
        worst = max(worst, abs(off[k - 1] ** 2 - cf2))
    M = mp.matrix(d, d)
    for k in range(d): M[k, k] = diag[k]
    for k in range(d - 1): M[k, k + 1] = M[k + 1, k] = off[k]
    ev = sorted(mp.re(x) for x in mp.eig(M)[0])
    cs = [mp.mpf(1)]
    for k in range(d): cs.append(cs[-1] * (k - d) * (k + U) / ((k + 1) * (k + V) * U))
    rts = sorted(mp.re(r) for r in mp.polyroots(cs[::-1], maxsteps=400, extraprec=400))
    err = max(abs(a - b) / V for a, b in zip(ev, rts))
    gers = max(abs(diag[k] - V) + (off[k - 1] if k > 0 else 0) + (off[k] if k < d - 1 else 0) for k in range(d)) / mp.sqrt(V * d)
    print("d=%d V=%d U=%d: max rel |eig - root| = %.1e, Gershgorin radius/sqrt(Vd) = %.4f, true max|root-V|/sqrt(Vd) = %.4f" % (
        d, V, U, err, gers, max(abs(r - V) for r in rts) / mp.sqrt(V * d)))
print("closed forms for diag - V and off^2: max abs error %.1e" % worst)
