# genus1777

Supplementary material for

> R. Auffarth, *A modular curve of genus 1777 with completely decomposable Jacobian*.

The article shows that the modular curve $X_\Gamma$ of the level-120 congruence subgroup

$$\Gamma=\{\gamma\in\mathrm{SL}_2(\mathbf Z):\ \gamma\equiv\lambda I \ (8),\ \ \gamma\equiv\pm\left(\begin{smallmatrix}1&*\\0&1\end{smallmatrix}\right)\ (3),\ \ \gamma\equiv\pm I\ (5)\}$$

has genus 1777 and a Jacobian isogenous over $\overline{\mathbf Q}$ to a product of elliptic curves.
The curve is isomorphic to $X_H(4800)$ with $H=\{d \equiv \pm1 \bmod 5\}$, so
$\mathrm{Jac}(X_\Gamma)\sim J_0(4800)\times J_{\chi_5}$, and the proof reduces to
newforms of trivial character and level dividing 4800 (Sections 2–6 of the article).

This repository contains everything needed to reproduce the computations of the proof,
their complete outputs, the data behind the tables, and the data of the independent
equivariant verification (§6.2).

## Contents

| Path | What it is | Paper |
|---|---|---|
| `scripts/verify_genus1777.py` | Main verification (PARI/GP via `cypari2`): checks (V0)–(V5) | App. B.1 |
| `scripts/crosscheck_modsym.py` | Second algorithm for computation (C1): modular symbols, checks (X1)–(X3) | App. B.1 |
| `scripts/verify_genus1777.m` | Magma version of computations (C1) and (C2), for readers with Magma (**not run by the author**) | App. B.1 |
| `scripts/make_catalogue.py` | Regenerates `data/newform_orbits.csv` | — |
| `outputs/verify_output.txt` | Complete output of `verify_genus1777.py` | App. B.1 |
| `outputs/crosscheck_output.txt` | Complete output of `crosscheck_modsym.py` | App. B.1 |
| `data/newform_orbits.csv` | All 396 Galois orbits of newforms of level M \| 4800, character 1 or χ₅: field, multiplicity σ₀(4800/M), traces of a_p for p < 100 | Table 2 |
| `data/table1_sixteen_orbits.csv` | The 16 non-rational trivial-character orbits, their CM or inner twist, Sturm level and bound, quaternion symbol | Table 1 |
| `data/table3_frobenius_traces.csv` | Frobenius traces of Jac(X_Γ) at 7 ≤ p ≤ 43 | Table 3 |
| `data/table4_newspaces.csv` | Dimensions and orbit counts of all nonzero newspaces S₂ⁿᵉʷ(M, ε) | Table 4 |
| `equivariant/` | Data and certificate of the independent equivariant verification | §6.2, App. B.2 |
| `SHA256SUMS` | Checksums of all files | — |

## What the checks establish

* **(V0)** All 70 newspace dimensions dim S₂ⁿᵉʷ(M, ε), M | 4800, ε ∈ {1, χ₅}, agree with an
  independent pure-Python evaluation of the Cohen–Oesterlé formula.
* **(V1)** dim S₂(4800, 1) = 889, dim S₂(4800, χ₅) = 888, total 1777.
* **(V2)** Every newspace is decomposed into Galois orbits, and the orbit degrees are checked to
  add up to its dimension, so the catalogue is complete. Exactly 16 trivial-character orbits are non-rational.
* **(V3)** For each of the 16, the CM or inner-twist identity is proved coefficient by coefficient
  through the Sturm bound (Lemma 3.5), in exact arithmetic.
* **(V4)** Consistency: every χ₅-orbit is a quartic twist of a trivial-character orbit.
* **(V5)** Consistency: the Frobenius traces of Table 3.
* **(X1)–(X3)** The orbit decomposition is recomputed with modular symbols, an independent algorithm.

## Requirements

* Python ≥ 3.10
* PARI/GP ≥ 2.17 with the `cypari2` bindings: `pip install -r requirements.txt`
  (the outputs here were produced with PARI/GP 2.17.2)
* About 8 GB of memory for the PARI stack

## How to run

```sh
cd scripts
python3 verify_genus1777.py  > ../outputs/my_verify_output.txt     # about 10 minutes
python3 crosscheck_modsym.py > ../outputs/my_crosscheck_output.txt  # about 15-20 minutes
magma verify_genus1777.m                                            # optional, Magma
```

Expected final lines:

```
All checks passed.
All modular-symbol cross-checks passed.
```

Each script stops with an assertion error if any check fails. Compare your output with the
files in `outputs/`.

## Checksums

```sh
sha256sum -c SHA256SUMS
```

## License

Scripts: MIT (see `LICENSE`). Data and outputs: CC BY 4.0 (see `LICENSE-DATA`).

## Citation

If you use this material, please cite the article (see `CITATION.cff`).
