# Equivariant verification (§6.2 of the article)

This folder holds the data of the independent verification sketched in §6.2. That method treats
X_Γ as a quotient of the full modular curve of level 360, decomposes its differentials into
44 determinant-twist blocks under GL₂(Z/360Z), bounds the conductors of the constituents with
Casselman's theorem, computes exact Hecke traces at the 354 primes 7 ≤ p ≤ 2399, and
identifies the constituents by exact rank certificates. It is logically independent of the
main proof in Sections 2–6.

## Files

* `blocks44.csv`: the 44 blocks. The columns are:
  * `block`: the block number.
  * `i8`, `i9`, `i5`: the representative, given as zero-based row indices in the exact character
    tables of GL₂(Z/mZ), m = 8, 9, 5, that ship with the certificate.
  * `d_C`: the block's contribution to the genus. These sum to 1777.
  * `g_rho`: the dimension of the representative's multiplicity space.
  * `candidate_columns_r`: the number of embedded candidate columns.
  * `minor_det_mod_1000000009`: the determinant of the recorded maximal minor modulo 10⁹ + 9.
    It is nonzero for every block.
  * `constituents_of_W_rho`: the newform orbits occurring in the representative, by LMFDB label,
    each with multiplicity one.

* `cert1777_full.json.gz`: the full certificate. Its SHA-256 digest is
  `7cb552e4d3ba8e03afcef1323df74e824e2ff80673e6bb9317e9f4669c4c42c6`, as stated in
  Appendix B.2 of the article. It contains:
  * the local character tables, the subgroup elements and the 44 blocks;
  * the candidate spaces and exact coefficient data;
  * the 15576 exact traces, the multiplicities and the maximal minors;
  * the coefficient fields and twist data of the non-rational forms.

* The Python verifiers for the certificate: group computations, trace formula, candidate lists,
  rank minors, twist and Brauer-class checks, and the property-first variant.

> **To do before publishing the repository:** add `cert1777_full.json.gz` and its Python
> verifiers to this folder, check that `sha256sum cert1777_full.json.gz` returns the digest above,
> and regenerate `../SHA256SUMS`.
