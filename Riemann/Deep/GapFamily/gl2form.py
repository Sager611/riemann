"""gapfamily: SDP form (with cuts) for the uniform measure at Lam, order J; writes logs/gl2_J<J>_<Lam>.json (cert5.py SRC format)."""
import sys, json
from sdpcut import sdp_cut
J = int(sys.argv[1]); lam = float(sys.argv[2].rstrip('m')); margin = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
m = [lam ** (2 * j) / (2 * j + 1) for j in range(J + 1)]
val, a, ph, mnv, rds = sdp_cut(m, J, rounds=14, margin=margin)
print('J=%d Lam=%s margin=%g: SDP value %.6g, cut rounds %d, min validity %.3g' % (J, sys.argv[2], margin, val, rds, mnv))
json.dump([{'lam': lam, 'val': val, 'minvalid': mnv, 'a%d' % J: [float(x) for x in a], 'phi%d' % J: {str(k): v for k, v in ph.items()}}],
          open('logs/gl2_J%d_%s.json' % (J, sys.argv[2]), 'w'), indent=1)
