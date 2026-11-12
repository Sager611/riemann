#!/usr/bin/env python3
"""gen2.py - generates Fail.lean (13 failing orders), Small.lean (exact sets for small N,
all 24 orders: witnesses and the sound search) and Main.lean (Theorem 1, Corollary).
Usage (repo root): python3 Riemann/Deep/Orders/Lean/gen2.py"""
import os, sys
sys.dont_write_bytecode = True  # keep ../__pycache__ out of the repo
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import verify as V

OPS = {'+': '.add', 'x': '.mul', '-': '.sub', '/': '.div'}
NAME = {'+': 'p', 'x': 'x', '-': 'm', '/': 'd'}
SYM = {'+': '+', 'x': '×', '-': '−', '/': '÷'}
nm = lambda o: ''.join(NAME[ch] for ch in o)
lo = lambda o: '(⟨%s⟩ : Order)' % ', '.join(OPS[ch] for ch in o)
sym = lambda o: ' '.join(SYM[ch] for ch in o)
spl = lambda l: '[' + ', '.join(map(str, l)) + ']'

def kfail(o):
    """k with no collection reachable in k unbounded moves (Proposition 1)."""
    if o[0] == '-': return 1
    if o[1] == '-': return 2
    assert o == '/x-+'; return 5

def seq(o, N):
    s = V.search(o, N); assert s is not None and V.check(o, s, N) is None
    return '[' + ', '.join('⟨%d, %d, %d⟩' % (x, y, v) for (x, op, y, v) in s) + ']'

def fail_file():
    out = ['import Riemann.Deep.Orders.Lean.Search', '',
           '/-! The 13 failing orders (NOTE.md §2, Proposition 1): no k-th move can make a new number. -/',
           'namespace Orders', '']
    for o in V.ORDERS:
        if V.CLAIM[o][1] is not None: continue
        k = kfail(o)
        out.append('/-- %s: no legal %d-move prefix, so a perfect N has N ≤ %d. -/' % (sym(o), k, k))
        out.append('theorem fail_%s : ∀ N ms, Perfect %s N ms → N ≤ 5 := fun N ms h =>' % (nm(o), lo(o)))
        out.append('  Nat.le_trans (bound_of_reach (k := %d) (by decide +kernel) N ms h) (by decide)' % k)
        out.append('')
    out.append('end Orders')
    return '\n'.join(out) + '\n'

def small_file(orders):
    out = ['import Riemann.Deep.Orders.Lean.Search', '',
           '/-! Exact sets of perfect N below N0 (working orders) or below 6 (failing orders):',
           'witnesses by the checker, non-existence by the sound search noPerfect. -/',
           'namespace Orders', '']
    for o in orders:
        spor, N0 = V.CLAIM[o]
        M = N0 if N0 is not None else 6
        out.append('theorem small_%s : ∀ N, N < %d → ((∃ ms, Perfect %s N ms) ↔ N ∈ %s) := by' % (nm(o), M, lo(o), spl(spor)))
        out.append('  intro N hN')
        out.append('  rcases (by omega : ' + ' ∨ '.join('N = %d' % N for N in range(M)) + ') with ' + ' | '.join(['rfl'] * M))
        for N in range(M):
            if N == 0:
                out.append('  · exact ⟨fun ⟨_, h⟩ => absurd h.1 (by omega), fun h => absurd h (by decide)⟩')
            elif N in spor:
                out.append('  · exact ⟨fun _ => by decide, fun _ => ⟨%s, (check_iff _ _ _).1 (by decide +kernel)⟩⟩' % seq(o, N))
            else:
                assert V.search(o, N) is None
                out.append('  · exact ⟨fun ⟨ms, h⟩ => absurd h (noPerfect_sound (by decide +kernel) ms), fun h => absurd h (by decide)⟩')
        out.append('')
    out.append('end Orders')
    return '\n'.join(out) + '\n'

def main_file():
    out = ['import Riemann.Deep.Orders.Lean.Fail'] + ['import Riemann.Deep.Orders.Lean.Small%d' % g for g in range(1, 7)]
    for o in V.DATA: out.append('import Riemann.Deep.Orders.Lean.D_%s' % nm(o))
    out += ['', '/-! Theorem 1 (exact sets of perfect N) and the Corollary (which orders work). -/',
            'namespace Orders', '']
    for o in V.ORDERS:
        spor, N0 = V.CLAIM[o]
        if N0 is not None:
            out.append('/-- Theorem 1, %s: perfect exactly for N ∈ %s or N ≥ %d. -/' % (sym(o), spl(spor), N0))
            out.append('theorem thm1_%s (N : Nat) : (∃ ms, Perfect %s N ms) ↔ N ∈ %s ∨ %d ≤ N := by' % (nm(o), lo(o), spl(spor), N0))
            out.append('  by_cases h : N < %d' % N0)
            out.append('  · rw [small_%s N h]; exact ⟨Or.inl, fun h\' => h\'.elim id (fun h\' => absurd h\' (by omega))⟩' % nm(o))
            out.append('  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_%s.works N (by omega)⟩' % nm(o))
        else:
            out.append('/-- Theorem 1, %s: perfect exactly for N ∈ %s. -/' % (sym(o), spl(spor)))
            out.append('theorem thm1_%s (N : Nat) : (∃ ms, Perfect %s N ms) ↔ N ∈ %s := by' % (nm(o), lo(o), spl(spor)))
            out.append('  by_cases h : N < 6')
            out.append('  · exact small_%s N h' % nm(o))
            out.append('  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_%s N ms hp) (by omega), fun h\' => ?_⟩' % nm(o))
            if spor:
                out.append('    simp only [List.mem_cons, List.not_mem_nil, or_false] at h\'; omega')
            else:
                out.append('    simp at h\'')
        out.append('')
    out += ['theorem cor_of_works {σ : Order} {R : Prop} (h : ∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms)',
            '    (hr : R) : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔ R := ⟨fun _ => hr, fun _ => h⟩', '',
            'theorem cor_of_fails {σ : Order} {R : Prop} (h : ∀ N ms, Perfect σ N ms → N ≤ 5) (hr : ¬R) :',
            '    (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔ R := by',
            '  constructor',
            '  · rintro ⟨N0, hN⟩',
            '    obtain ⟨ms, hms⟩ := hN (N0 + 6) (by omega)',
            '    have := h _ _ hms; omega',
            '  · intro r; exact absurd r hr', '']
    for o in V.ORDERS:
        spor, N0 = V.CLAIM[o]
        R = '(%s.s3 = .sub ∨ %s.s4 = .sub) ∧ %s ≠ ⟨.div, .mul, .sub, .add⟩' % (lo(o), lo(o), lo(o))
        out.append('theorem cor_%s : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect %s N ms) ↔ (%s) :=' % (nm(o), lo(o), R))
        if N0 is not None:
            out.append('  cor_of_works ⟨%d, D_%s.works⟩ (by decide)' % (N0, nm(o)))
        else:
            out.append('  cor_of_fails fail_%s (by decide)' % nm(o))
        out.append('')
    out += ['/-- **Corollary.** An order σ of + × − ÷ has perfect sequences for all large N iff − is',
            'third or fourth and σ ≠ ÷ × − + (11 of the 24 orders). -/',
            'theorem corollary (σ : Order) (hσ : σ.IsPerm) :',
            '    (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔',
            '      ((σ.s3 = .sub ∨ σ.s4 = .sub) ∧ σ ≠ ⟨.div, .mul, .sub, .add⟩) := by',
            '  obtain ⟨a, b, c, d⟩ := σ',
            '  cases a <;> cases b <;> cases c <;> cases d <;>',
            '  first',
            '  | exact absurd hσ (by unfold Order.IsPerm; decide)']
    for o in V.ORDERS:
        out.append('  | exact cor_%s' % nm(o))
    out += ['', 'end Orders']
    return '\n'.join(out) + '\n'

if __name__ == '__main__':
    jobs = [('Fail', fail_file), ('Main', main_file)]
    for g in range(6):
        jobs.append(('Small%d' % (g + 1), (lambda g=g: small_file(V.ORDERS[4 * g:4 * g + 4]))))
    for fn, f in jobs:
        txt = f(); open(os.path.join(HERE, fn + '.lean'), 'w').write(txt)
        print(fn, len(txt.splitlines()), 'lines', len(txt), 'bytes')

