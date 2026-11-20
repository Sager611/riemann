# Reviewer verdict (one independent reviewer subagent, verdict only, evidence pasted inline; 20 Nov 2026)

Brief: claims C1–C7 (target footprint; host-step failure; Lemma A; Lemma B; rewrites R1–R3; transferring steps; NOGO wording) with the Lean quotes and probe numbers inline (.scratch/alpext/ not kept; the claims are those of SCOPE.md §0–§2 and PROBES.md).

Verdict (verbatim, abridged):
- C1 CANNOT JUDGE: the quoted interfaces establish target-entry points, not exhaustiveness; absence of PChar does not exclude indirect dependencies.
- C2 AGREE, for the stated host mechanism: a cyclotomic host cannot realize irreducible ρ; realizing ψ faithfully over F forces H ⊂ E, impossible for CM E. Literally, PChar rejects ψ earlier.
- C3 AGREE (Lemma A; K of signature (1,1) excludes K and H from every CM field).
- C4 AGREE (Lemma B), reading L_ψ as the field cut out by ψ; ring-class conjugation acts by inversion.
- C5 DISAGREE as stated: the CM and fixed-host orthogonality obstructions hold; centrality is necessary uniformly over class weights, not for every particular weight; unit-valued twists do not exclude new trace-sensitive constructions.
- C6 DISAGREE as stated: finite checks establish no uniform divisibility theorem or wholesale transfer; coefficient positivity requires −1/2 ≤ w ≤ 1 unless w ≥ 0 is assumed.
- C7 DISAGREE with the wording: ψ's values lie in the CM field ℚ(ζ₃); say that the field H cut out by ψ embeds in no CM field.
- Overall: 'BLOCKED at the host step': AGREE for the unchanged archive mechanism; not a proof that every twisted or representation-valued replacement fails. Full package INCOMPLETE. Most serious gap: a uniform inert-prime divisibility statement (C6) is only checked finitely.

Fixes made:
- C1: SCOPE §0 now says the core depends on ψ only through the host it is handed (indirect dependence = the host step itself), at the type level; exhaustiveness beyond types is listed under Limits.
- C5: the partner argument is stated for all class weights (AskMainBalance quantifies ∀ w, FirstPrimeNeSecondPrime.lean:99); R3 is marked FAILS for the archive's mechanism and NEEDS-NEW-INPUT (UNKNOWN) for a new trace-sensitive construction.
- C6: P2 is labelled a finite check (classical only D^p ≡ 0 mod P, Hasse invariant 0); positivity stated for 0 < w ≤ 1 (CharFamily weights are positive, CharFamily.lean:31).
- C7: NOGO line reworded ('the field cut out by ψ').
- Verdict line: BLOCKED for the archive's route; open (UNKNOWN) for mechanisms not in the archive.
