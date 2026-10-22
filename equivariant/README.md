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

* `repro1777/`: the full certificate package, unchanged from the author's reproducibility bundle
  (see `repro1777/README.md` for its own description). It contains:
  * `repro1777/cert1777_full.json.gz`: the certificate. Its SHA-256 digest is
    `7cb552e4d3ba8e03afcef1323df74e824e2ff80673e6bb9317e9f4669c4c42c6`, as stated in Appendix B.2
    of the article. It contains the local character tables, the subgroup elements and the 44 blocks;
    the candidate spaces and exact coefficient data; the 15576 exact traces, the multiplicities and
    the maximal minors; and the coefficient fields and twist data of the non-rational forms.
  * Three verifiers that need only the Python standard library and read only the certificate:
    ```
    cd repro1777
    python3 verify1777.py      cert1777_full.json.gz   # V1-V5   (a few minutes)
    python3 verify_traces.py   cert1777_full.json.gz   # V6      (a few minutes)
    python3 verify_complete.py cert1777_full.json.gz   # V7a, V7b
    ```
  * The scripts `s1_tables.g` to `s10_export.py` and the GAP tables `tables_*.txt`, which regenerate
    the certificate from scratch (GAP 4, PARI/GP, Python with numpy and sympy).
* `verifier_output.txt`: the output of the three verifiers above. All checks pass.
