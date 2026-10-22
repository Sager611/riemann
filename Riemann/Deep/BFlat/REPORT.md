NO DISCOVERY: B♭ is pinned by one tight inequality, the endpoint of the critical row class (δ = 0.4693, x = 1/2): its count 0.590332 is set jointly by the Möbius second moment (prime cost 2), the plain fourth moment (prime cost 6κ) and sixth-power amplification (5/6), and its numerator by the fourth moment; c/12, 2/27, 1/13, Lemmas 18.2/18.3 and the 13/15 floor are all slack.

# Lowering the quasi-RH constant B♭ = 0.874886647: what binds (2026-10-22)
Details: [NOTE.md](NOTE.md). Scripts: bflat_exact.py, bflat_sens.py, bflat_checks.py (each < 7 s, < 100 MB).

## Records
- Nothing below B♭ on arXiv (math.NT October listing, new submissions, API), Zenodo (32 quasi-RH and 31 half-plane records), Hexagon (all 342 works) or openai/math (last commit 12 Oct). The lowest other value is Newton's 0.87488664962 (Zenodo), which Liu cites.

## Liu's stage map
- **One step.** A proved bound β_* ≤ B licenses κ = 2B − 1 and yields every B′ > F♭(κ) = 11/12 − ℓ*(κ)/4. Here ℓ* comes from the critical class δ_c(κ), the root of cr(δ,κ) = 0, where the endpoint does not depend on b. Full admissibility list: NOTE §1.
- **Fixed point.** The root of Liu's quartic 133488B⁴ − 376596B³ + 384534B² − 167541B + 26231 (resultant reproduced exactly): B♭ = 0.8748866470399641. The map contracts by 0.003; iterates from 7/8 miss B♭ by 3.4e-7, 1.0e-9, 3.1e-12, ….
- **κ barely matters:** F♭(1) = 0.875169, F♭(3/4) = 0.874887, F♭(1/2) = 0.874312 (formal).

## What binds
- **One tight inequality:** E(δ_c, 1/2) = 0, a double root at δ_c = 0.469345, with ℓ* = 0.167120, b* = 0.086535.
- **Count there:** R = 0.590332 (density hypothesis: 0.530655). Three estimates are tight at once: the Möbius moment (length 0.746, primes 0.127), the plain fourth moment (length 0.418, primes 0.036) and amplification (cutoff 1.164). The numerator is (3R+1)/4 with shift g/(6κ); the amplitude cap x ≤ 1/2 is active.
- **Sensitivity of B♭ per unit change:** amplification exponent 5/6: +0.0074; Möbius prime cost 2: +0.0021; plain prime cost 6: +0.00019; numerator moment order 4: −0.0035.
- **Slack at B♭:** c/12 (0.036 vs 0.074), 1/13 (0.0741 vs 0.0769), prime supply (0.177 vs 0.210), 13/15 floor (ℓ ≤ 1/5, slack 0.033), Lemmas 18.2/18.3 (needed, but no constant), 2/27 (OpenAI's bound for κ ≥ 3/4, unused here).

## Attempts, none of which beats B♭
- **More stages:** the map contracts to B♭. **Numerator chosen class by class:** at ℓ* + ε the worst class stays at +1.089ε.
- **Recombining the three moment lemmas** (Hölder, several cutoffs, a trivial second factor, primes alone, a Heath-Brown split giving 0.666, amplification with primes, the general large sieve, Halász): all worse.
- **Constants inside the lemmas:** 5/6 and 6κ come from the order 6 of the characters (every sixth-power multiplier is already used; the budget M/6 is what the 5M/6 threshold leaves). The Möbius cost 2 comes from the mark invariants of OpenAI §17, which never bound a prime sum pointwise, so β_* cannot enter.

## What a record needs
Any one new theorem for the thin sextic family; none is known: a count below U^0.590332 at (0.4693, 1/2) (a density exponent below 2.2249 at σ = 0.7347), amplification below 5/6, a Möbius prime cost below 2, a plain prime cost below 6κ, or a numerator average below U^{(3R+1)/4−θ/2}.

Suggested NOGO line: "B♭ 0.874886647: one critical class binds (count via M1 cost 2, M2 cost 6κ, amplification 5/6; numerator 4th moment); c/12, 2/27, 1/13, L18.2/18.3, 13/15 slack; needs a new sextic-family estimate."
