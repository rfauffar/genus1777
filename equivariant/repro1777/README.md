# Genus-1777 completely decomposable modular curve: reproducible certificate

## Statement

Let

Γ = { γ ∈ SL₂(Z) : γ ≡ λI (mod 8) with λ² ≡ 1,  γ ≡ ±(1 b; 0 1) (mod 3),  γ ≡ ±I (mod 5) }.

Then X_Γ has genus 1777, and Jac(X_Γ) is isogenous over Q̄ to a product of elliptic curves.

Γ has level 120, and contains Γ(120) and −I. Its image Γ̄ in SL₂(Z/360) is G₈ × G₉ × G₅, of orders 4 × 162 × 2. The mod-3 condition uses upper-unipotent matrices; this is conjugate to the lower-unipotent version used in earlier files.

## Verifying (no GAP or PARI needed)

The three verifiers need Python ≥ 3.8 and only its standard library. They read only `cert1777_full.json.gz`, which also contains the GAP character-table output.

```
python3 verify1777.py     cert1777_full.json.gz     # V1-V5    (~2 min)
python3 verify_traces.py  cert1777_full.json.gz     # V6       (~4 min)
python3 verify_complete.py cert1777_full.json.gz    # V7a, V7b (~1 min)
```

| Check | What is verified |
|---|---|
| V1 | Coset enumeration of Γ̄ in SL₂(Z/360) from the stored generators: index 23040, e₂ = e₃ = 0, 288 cusps, genus 1777 |
| V2 | The 44 block dimensions sum to 1777 |
| V3 | For every block ρ and all 354 primes 7 ≤ p ≤ 2400, p ∤ 360: Σ_f n_f · Tr_{K_f/Q(χ)} a_p(f) = −N_p(ρ), exactly in Q(ζ₂₄) |
| V4 | For every block, the matrix of candidate coordinate vectors at the recorded witness rows has nonzero determinant mod ℓ = 1000000009, with ζ₂₄ ↦ the recorded root of unity. So the candidate columns are linearly independent, and the multiplicities in V3 are the unique decomposition |
| V5 | Every non-rational newform orbit used is elliptic over Q̄ (details below) |
| V6 | The traces −N_p(ρ) are recomputed from the GAP character tables, with class numbers counted here from reduced binary quadratic forms, and agree with the certificate |
| V7a | The conductor bounds and central character of each block representative are recomputed from the tables. The candidate spaces are exactly all (M, ε) they allow: M = 2^a·3^b·5^c within the bounds, M ≥ 11, ε ∈ {ω, ω̄} factoring through (Z/M)^\* |
| V7b | For each of the 38 candidate spaces, the newspace dimension from the Cohen–Oesterlé formula (with the new-subspace inclusion–exclusion β(p) = −2, β(p²) = 1) equals both the recorded dimension and the number of newforms listed |

V5 covers three kinds of orbit:
- **CM orbits** (D = −4, −20): a_n(1 − χ_D(n)) = 0 for all n ≤ Sturm(lcm(M, D²)) with gcd(n, MD) = 1.
- **Inner-twist orbits:**
  - each recorded automorphism σ(y) is a root of the Hecke field polynomial, and there are [K_f : Q] of them;
  - σ(a_n) = χ_D(n)·a_n holds for n up to the Sturm bound;
  - Quer's square classes of a_p²/ε(p) and the matching twist discriminants are recomputed from the q-expansions, with ε taken from the stored character;
  - all Hilbert symbols are recomputed, and every Brauer class is split.
- **Rational orbits** need no check.

## Regenerating everything (GAP 4.x, PARI/GP ≥ 2.15 via cypari2, Python 3 with numpy and sympy)

```
gap -q s1_tables.g          # character tables of GL2(Z/8), GL2(Z/9), GL2(Z/5) -> tables_*.txt
python3 s2_classes.py       # conjugacy classes matched to GAP (sizes, orthogonality checked) -> classes.pkl, tables.pkl
python3 s3_group.py         # Gamma-bar, coset genus, exact Chevalley-Weil multiplicities -> group.pkl
python3 s4_blocks.py        # 44 determinant-twist blocks -> blocks.pkl
python3 s5_bounds.py        # local conductor bounds (corrected lemma) -> bounds.pkl
python3 s6_traces.py        # exact traces in Q(zeta_24), 354 primes -> traces.pkl
python3 s7_newforms.py      # candidate spaces, exact Hecke data from PARI -> spaces/*.pkl, plan.pkl   (rerun until done)
python3 s8_match.py         # rank witnesses, multiplicities, exact identities -> match.pkl
python3 s9_elliptic.py      # CM / inner-twist / Brauer certificates -> elliptic.pkl   (needs used_orbits.pkl, see s8)
python3 s10_export.py       # writes cert1777_full.json.gz
```

Between s8 and s9, `used_orbits.pkl` comes from the multiplicities in `match.pkl`; s8's logs show how. The exact arithmetic in s3–s6 uses Z[ζ₂₄] modulo Φ₂₄ (`cyclo.py`).

## What is assumed

The verifiers check the arithmetic. The conclusion also rests on the following standard results.

1. **Decomposition into blocks.** J(X_Γ) ~ ⊕_ρ V_ρ^{dim ρ^Γ̄} over Q̄, where ρ runs over irreducible representations of GL₂(Z/360). This is the adelic Eichler–Shimura relation for X(360) together with Chevalley–Weil. Here dim V_ρ = g_ρ, computed by the Chevalley–Weil formula with the elliptic and cuspidal stabilizers (s3).
2. **Traces.** The trace of Frobenius at p on V_ρ is −N_p(ρ), by the Lefschetz/Eichler–Selberg trace formula in the class-number form implemented in s6 and V6.
3. **Newforms in each block.** V_ρ is a sum of pieces A_f indexed by newforms f whose local representations π_{f,p} contain ρ_p. By Casselman's newform theory and the conductor lemma below, a(π_{f,p}) is bounded and the nebentypus is ω or ω̄. Hence the candidate spaces of V7a are complete.
4. **Uniqueness.** With linearly independent candidate columns (V4), the exact identity (V3) determines the multiplicities.
5. **Ellipticity.**
   - Rational newforms give elliptic curves over Q (Eichler–Shimura).
   - CM newforms give A_f ~ E^d over Q̄ (Shimura).
   - Non-CM newforms whose inner twists exhaust Gal(K_f/Q) have F = Q (Ribet). Then A_f ~ E^d over Q̄ if and only if the Brauer class is split, and that class is computed by Quer's formula.
6. **Sturm bounds.** Twist and CM identities are checked for n ≤ Sturm(lcm(M, D²)) with gcd(n, MD) = 1. This is the standard convention for comparing a newform with its twist. A version that zeroes out the bad coefficients first would need level (MD)², and has not been done.

**Conductor lemma**, as corrected after review. Let ρ be an irreducible representation of GL₂(Z/p^e), and let n ≥ 1, 0 ≤ k ≤ e, 0 ≤ n − k ≤ e. Put

H_{k,n} = { (a, p^k b; p^{n−k} c, d) mod p^e : det a unit, d ≡ 1 mod p^{min(n,e)} }.

If ρ^{H_{k,n}} ≠ 0, then every irreducible infinite-dimensional π with ρ ⊂ π^{K(p^e)} has a(π) ≤ n.

*Proof.* For δ = diag(p^k, 1) we have δK₁(pⁿ)δ⁻¹ ⊆ the preimage of H_{k,n}. Here n ≥ 1 forces c′ ≡ 0 mod p, so a′ is a unit. Hence π(δ⁻¹)v is K₁(pⁿ)-fixed, and Casselman gives a(π) ≤ n. The case n = 0 is the spherical case: H = GL₂(Z/p^e), and only the trivial representation has invariants.

## What is trusted from software, and how it is cross-checked

- **GAP character tables.** The verifiers check that the GAP classes cover the group, that class sizes agree with recomputed conjugacy orbits, and (in s2) orthogonality.
- **PARI eigenforms** (the q-expansion coefficients). The verifiers check that they satisfy the exact trace identities at all 354 primes (V3), that the newspace dimensions match Cohen–Oesterlé (V7b), and the twist and CM identities (V5).
- What the verifiers do **not** do is recompute the eigenforms themselves. An independent check would compare the stored a_p of the 49 orbits with Magma, Sage or the LMFDB. The orbits have levels 15 to 2400, with nebentypus trivial or quadratic.

## Changes from the earlier certificate

The corrected conductor lemma gives smaller candidate sets. The decomposition is now expressed with 49 newform orbits of level at most 2400:
- 36 rational;
- 2 CM;
- 11 with inner twists and F = Q, all with split Brauer class.

The earlier certificate used 61 orbits at levels up to 4800. Both describe the same abelian variety through different block representatives (twists).
