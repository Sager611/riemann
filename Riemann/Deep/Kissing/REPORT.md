RECORD: kissing number in dimension 38 ≥ 624020 (was 611522)

- New configuration: Kravatskiy's 611522-point configuration in R^38 (current record, Cohn's table 4 Nov 2026) with its
  equator (Lin's 6855 lines = 13710 points) replaced by 13104 lines = 26208 points. Gain +12498.
- The 13104 lines: all norm-48 vectors of M = {x in Leech : F x = 0 mod 3}, a sublattice of index 3^6 with no
  norm-32 vector, so any two differ and sum to vectors of norm >= 48, i.e. |v.w| <= 24 (cos <= 1/2).
- Found via an Eisenstein structure of the Leech lattice: L/theta L = F_3^12, the minimal-vector classes have a
  two-valued Fourier transform (-90, 639), which forces dim <= 6 for avoiding subspaces; a greedy search finds 6.
- Lin's local-search value was 6855 (alpha of his 8.4M-vertex graph known only in [6855, 98280]); ours 13104;
  it is locally rigid (every other line meets >= 11 of ours); LP bound 40320 lines.
- Verified: verify.py (stdlib, exact, 3.6 min, includes brute force over all 85.85M pairs); Lin's own C checker
  (pairs and all 196560 minimal vectors); Kravatskiy's two verifiers with only the equator swapped (both pass 624020).
- Task 1: re-verified Kravatskiy's 611522 (his verifiers), Lin's 6855 equator, and Lin's 370240 points in R^33 by
  brute force over all 6.85e10 pairs (max inner product 16 = |x|^2/2).
- Both target papers are already superseded (Kravatskiy repo: 32-39 per Cohn's table; 49-60 with a 7783-line P48p
  class). 49D would need a 7784-line class; not reached.
- Caveat: non-equator points are Kravatskiy's certified data (commit b814435); only the equator is new.
- Files: Riemann/Deep/Kissing/ (NOTE.md details, verify.py, equator_lines_13104.txt, F3.txt, golay_basis.txt,
  omega8.txt, discover.py, construct.py, replay_kravatskiy.py).

