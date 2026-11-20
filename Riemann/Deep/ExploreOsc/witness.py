# witness.py (explore60): deterministic members of the class S (flat spectral measure on [-pi, pi], unfolded units: maximal
# zero density 1 per unit) with few zeros. Amplitude-modulated chirp: instantaneous frequency w(u) spends time density
# tau(w) ~ 1/(w + delta) at w, amplitude a(w) = sqrt(w + delta) so that a^2 tau is flat; zeros = zeros of cos(phase).
# Up-sweep then mirrored down-sweep, repeated. Prints the power fraction below a*pi (flat: a) and zeros per unit length.
import sys, numpy as np
delta = float(sys.argv[1]); K = float(sys.argv[2]); reps = int(sys.argv[3]); du = 0.02
U = K*np.log((np.pi + delta)/delta)                 # one up-sweep: w(u) = delta (e^{u/K} - 1), w(U) = pi
u = np.arange(0.0, U, du)
w = delta*(np.exp(u/K) - 1)
ph = delta*K*(np.exp(u/K) - 1) - delta*u            # phase = int_0^u w
amp = np.sqrt(w + delta)
def piece(sign):                                    # up-sweep followed by its mirror image (continuous phase)
    up = amp*np.cos(ph); dn = (amp*np.cos(2*ph[-1] - ph))[::-1]
    return np.concatenate([up, dn])
x = np.concatenate([piece(1) for _ in range(reps)])
P = np.abs(np.fft.rfft(x*np.hanning(len(x))))**2
f = 2*np.pi*np.fft.rfftfreq(len(x), du)
fr = [P[f < a*np.pi].sum()/P[f < 1.5*np.pi].sum() for a in (0.1, 0.25, 0.5, 0.75, 1.0)]
z = int((np.sign(x[:-1])*np.sign(x[1:]) < 0).sum())
L = len(x)*du
pred = 1/np.log((np.pi + delta)/delta) - delta/np.pi
print('AM chirp delta=%g K=%g period=%.1f length=%.0f: F(.1,.25,.5,.75,1)=%s  zeros/length=%.4f (predicted %.4f; max 1)' %
      (delta, K, 2*U, L, ' '.join('%.3f' % v for v in fr), z/L, pred))
# windowed flatness at scale 4U
Wn = int(4*U/du); out = []
for i in range(0, len(x) - Wn, Wn):
    seg = x[i:i + Wn]*np.hanning(Wn); Ps = np.abs(np.fft.rfft(seg, 4*Wn))**2; fs = 2*np.pi*np.fft.rfftfreq(4*Wn, du)
    out.append([Ps[fs < a*np.pi].sum()/Ps[fs < 1.5*np.pi].sum() for a in (0.25, 0.5, 0.75)])
out = np.array(out)
print('  windows of length %.0f (n=%d): F(.25,.5,.75) min/max = %s' % (4*U, len(out),
      ' '.join('%.3f/%.3f' % (out[:, j].min(), out[:, j].max()) for j in range(3))))
