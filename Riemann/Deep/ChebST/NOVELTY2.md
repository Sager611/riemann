NEW: no external source states, for a fixed non-CM newform f, primes in every ideal class of Q(√D) with cos θ_p in an interval J (or a Sato–Tate discrepancy bound inside classes, or sectors) with any power of ∣D∣ — for almost all D, on average over D, or per D unconditionally. Only weaker relatives are known (table). The only prior statement of this shape is the project's own RSFam entry (signs, degree ≤ 2; D_M sketched).

# ChebST novelty audit 2 (adversarial, 2026-10-21)
Scratch `.scratch/stnov/` (deleted after the run). Goal: find what [NOVELTY.md](NOVELTY.md) missed.

## Searched (new ground)
- arXiv math.NT monthly listings 2025-01 … 2026-10: 22 pages, 8,298 entries including cross-lists. Abstracts of all 8,298 fetched through the API (id_list) and scanned locally (Sato–Tate / Hecke eigenvalue / Frobenius × class group / ideal class / quadratic form / Chebotarev / zero density / large sieve): 26 Sato–Tate abstracts, 49 combined hits, none relevant.
- arXiv API: 43 author queries (Thorner, Zaman, Lemke Oliver, Pierce, Turnage-Butterbaugh, M. M. Wood, P.-J. Wong, M. R. and V. K. Murty, E. Chen, P. Park, Swaminathan, Humphries, Pascadi, Kowalski, Michel, Bilu, Deshouillers, Gun, Naik, Akbary, J. Park, Zywina, Brumley, Cureton). Also 48 new phrase queries (split primes, norm form, Hilbert class field, Satake × class group, Rankin–Selberg × class-group/ring-class characters, Gaussian primes, every/each ideal class, ray class, Siegel zero, exceptional set, theta lift, function fields, …) and every "Sato Tate" abstract of all years (170). Nothing found.
- Forward citations: OpenAlex for 22 seed papers (Lemke Oliver–Thorner, Thorner 2021/2025, Wong 2019, Thorner–Zaman, Zaman, BTZ, Humphries–Thorner, PTBW, Kowalski–Michel, Duke–Kowalski, Bucur–Kedlaya, Bucur–Fité–Kedlaya, V. K. Murty, Harris, Duke–Friedlander–Iwaniec, Malik–Prabhu, …): 439 citing works with abstracts. Semantic Scholar citations of 11 arXiv papers: 210 works. 32 keyword hits in total, all already known.
- Zenodo: 11 boolean queries, 330 records (all 87 "Sato-Tate"). Only progressions (Lee), Caldin, R. Chen.
- openai/math: 749 preprints. GitHub code search of their TeX/Lean: "Sato-Tate" and "Sato Tate" give 0 hits. GitHub-wide search for "Sato-Tate" with "ideal class" / "class group" / "binary quadratic form": 0 hits.
- Hexagon: all 342 works (2609.00001–2610.00441), titles and abstracts: no Sato–Tate result.
- Full text (TeX) read: Frei–Sofos–Widmer 2610.18491, Bucur–Kedlaya 1301.0139, BTZ 1804.06402 (no Sato–Tate, no ideal classes). GXL evidence API: HTTP 429 again.

## What is known (weaker than the claim)
| Statement | Status | Source |
|---|---|---|
| fixed D, qualitative: θ_p ST-distributed over the primes of one class | known: special case of joint Chebotarev–Sato–Tate for the ring class field | Wong, JNT 196 (2019) (from BLGHT) |
| every D, polynomial in ∣D∣, under GRH | routine but never written for classes (our estimate: ∣D∣^{1+ε} for fixed J) | effective Chebotarev–Sato–Tate: V. K. Murty 1985, Bucur–Kedlaya 1301.0139, Wong 2019 |
| every D, unconditional, huge inexplicit power | derivable, stated nowhere (RSFam caveat); BTZ's density is now unconditional | BTZ 1804.06402 + Pascadi–Thorner 2508.14888 item (3) |
| no Sato–Tate condition (m = 0), every D | known | Thorner–Zaman 1604.01750 (∣D∣^694); Gaudet 2408.01610 (D prime) |
| progressions instead of classes, almost all moduli | analogue known | Bombieri–Vinogradov for GL_n (Jiang–Lü–Thorner–Wang 2104.02711) |
| signs, degree ≤ 2 minorants, almost all D | ours, not external | RSFam (DISCOVERIES) |

## New context missed before: 2026 quasi-RH claims
- These claims are unrefereed and backed by Lean. They give zero-free half-planes for Dirichlet L-functions (OpenAI, 7/8), for all finite-order Hecke L-functions of every number field (Cureton–Thorner, Hexagon 2610.00200: Re s > 31/32) and for all unitary Weil-group L-functions, Grössencharacters included (Hexagon 2610.00216: 7/8). Frei–Sofos–Widmer 2610.18491 apply them, conditionally, to Chebotarev.
- If these claims hold, the case m = 0 (and sectors without a Sato–Tate condition, and CM forms) becomes polynomial in ∣D∣ for every D. None of them covers L(Sym^m f ⊗ θ_χ) for non-CM f: Hexagon 2610.00178 says explicitly that no automorphic statement is inferred. So Theorems 1–3 are untouched, but the NOTE should say in one sentence that its m = 0 input is superseded per D if these claims hold.

## Verdict and residual risk
NEW for Theorems 1–3 as stated: almost all D, explicit A(J), discrepancy, sectors. It is an increment over RSFam, so merging the two entries (REPORT's suggestion) is reasonable.
Residual risk: arXiv totals look low (170 "Sato Tate" abstracts across all years), so API coverage may be incomplete. Semantic Scholar keyword search was rate-limited (1 of 7 queries ran). GXL was unavailable. Wong 2019 was seen at abstract level only. Journals outside arXiv were seen only through OpenAlex/S2 metadata.
