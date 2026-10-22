#!/usr/bin/env python3
"""
Independent verification for
  "A modular curve of genus 1777 with completely decomposable Jacobian".

Requires PARI/GP through the cypari2 bindings (tested with PARI/GP 2.17.2):
    pip install cypari2
Runtime: a few minutes on a laptop.

Checks performed (numbering as in Appendix B of the paper):
  (V0) Every newspace dimension dim S_2^new(M, eps), M | 4800, eps in {1, chi_5}, computed
       by PARI agrees with an independent pure-Python evaluation of the Cohen-Oesterle
       dimension formula followed by Moebius-type inversion (no PARI involved).
  (V1) dim S_2(Gamma_0(4800)) + dim S_2(4800, chi_5) = 1777, and the same total is
       recovered from the newform decomposition with multiplicities sigma_0(4800/M).
  (V2) For every M | 4800 and eps in {1, chi_5}: the newspace is decomposed into Galois
       orbits, and the sum of the coefficient-field degrees of the orbits is checked to
       equal the dimension of the space (so the list of orbits is complete).  Lists the
       trivial-character orbits with K_f != Q.
  (V3) For each of the 16 trivial-character orbits with K_f != Q: a CM self-twist or an
       inner twist (sigma, chi_D), proved by comparing coefficients through the Sturm
       bound of level L = lcm(M, D^2), coefficient by coefficient in exact arithmetic;
       the Hilbert symbols (t, D)_v are evaluated exactly at v = oo, 2, 3, 5 (the only
       places where they can be nontrivial) as a redundant check of the norm identities.
  (V4) (Consistency only.)  Every chi_5 newform twisted by psi^{-1}, psi of order 4
       mod 5, matches a trivial-character newform of level dividing 4800 at 7 <= p < 400.
  (V5) (Consistency only.)  Frobenius traces of the Jacobian at 7 <= p <= 43.
"""
import math
from fractions import Fraction
import cypari2

pari = cypari2.Pari()
pari.allocatemem(8 * 10**9)

N = 4800
DIVS = [d for d in range(1, N + 1) if N % d == 0]
QUAD = [-3, -4, 5, -8, 8, 12, -15, -20, 24, -24, 40, -40, 60, 120, -120]


def sigma0(n):
    return sum(1 for d in range(1, n + 1) if n % d == 0)


def kron(D, n):
    return int(pari.kronecker(D, n))


def sturm_bound(L):
    """ceil( (k/12) [SL2(Z):Gamma_0(L)] ) for k = 2."""
    num, den = L, 6
    for p in range(2, L + 1):
        if L % p == 0 and pari.isprime(p):
            num *= p + 1
            den *= p
    return -(-num // den)


# ---------------------------------------------------------------------------
# Independent dimension formula (pure Python, no PARI): Cohen-Oesterle, weight 2.


def _factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def _eps(chi, x):
    return 1 if chi == 1 else {0: 0, 1: 1, 4: 1, 2: -1, 3: -1}[x % 5]


def co_dim_S2(n, chi):
    """dim S_2(n, eps) for eps = 1 or chi_5 (5 | n), by the Cohen-Oesterle formula."""
    f = _factor(n)
    psi = Fraction(n)
    for p in f:
        psi *= Fraction(p + 1, p)
    lam = 1
    for p, r in f.items():
        s = 1 if (chi == 5 and p == 5) else 0
        if 2 * s <= r:
            lam *= (p ** (r // 2) + p ** (r // 2 - 1)) if r % 2 == 0 else 2 * p ** (r // 2)
        else:
            lam *= 2 * p ** (r - s)
    s4 = sum(_eps(chi, x) for x in range(n) if math.gcd(x, n) == 1 and (x * x + 1) % n == 0)
    s3 = sum(_eps(chi, x) for x in range(n) if math.gcd(x, n) == 1 and (x * x + x + 1) % n == 0)
    val = psi / 12 - Fraction(lam, 2) - Fraction(s4, 4) - Fraction(s3, 3) + (1 if chi == 1 else 0)
    assert val.denominator == 1
    return int(val)


def co_dim_new(n, chi):
    """Newspace dimension: sum over M | n of beta(n/M) dim S_2(M), beta = mu * mu."""
    def beta(m):
        out = 1
        for _, k in _factor(m).items():
            out *= {1: -2, 2: 1}.get(k, 0)
        return out
    return sum(beta(n // M) * co_dim_S2(M, chi) for M in range(1, n + 1)
               if n % M == 0 and (chi == 1 or M % 5 == 0))


def spaces():
    for M in DIVS:
        for chi in ([1, 5] if M % 5 == 0 else [1]):
            mf = pari(f"mfinit([{M},2,{chi}],0)")
            if int(pari.mfdim(mf)) > 0:
                yield M, chi, mf


def main():
    # (V0)
    nsp = 0
    for M in DIVS:
        for chi in ([1, 5] if M % 5 == 0 else [1]):
            nsp += 1
            a = co_dim_new(M, chi)
            b = int(pari.mfdim(pari(f"mfinit([{M},2,{chi}],0)")))
            assert a == b, (M, chi, a, b)
    print(f"(V0) {nsp} newspace dimensions agree with the Cohen-Oesterle formula")

    # (V1)
    g1 = int(pari("mfdim([4800,2,1],1)"))
    g5 = int(pari("mfdim([4800,2,5],1)"))
    print(f"(V1) dim S2(G0(4800)) = {g1}, dim S2(4800,chi5) = {g5}, sum = {g1 + g5}")
    assert g1 + g5 == 1777

    # (V2)
    total, nonrat_triv = 0, []
    counts = {}
    for M, chi, mf in spaces():
        m = sigma0(N // M)
        degs = [int(pari.poldegree(P)) for P in pari.mffields(mf)]
        assert sum(degs) == int(pari.mfdim(mf)), (M, chi)   # the orbit list is complete
        print(f"     S2new({M},{'1' if chi == 1 else 'chi5'}): dim {int(pari.mfdim(mf))}, "
              f"orbit degrees {sorted(degs)}, multiplicity {m}")
        for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
            d = int(pari.poldegree(P))
            total += m * d
            key = (chi, str(pari.polredabs(P)) if d > 1 else "Q")
            counts[key] = counts.get(key, 0) + 1
            if chi == 1 and d > 1:
                nonrat_triv.append((M, F, P))
    print(f"(V2) weighted total = {total}; orbit counts by (character, field): {counts}")
    assert total == 1777 and len(nonrat_triv) == 16

    # (V3)
    for M, F, P in nonrat_triv:
        a0 = pari.mfcoefs(F, 400)
        good = [p for p in range(2, 400) if pari.isprime(p) and M % p]
        s = [x for x in pari.nfgaloisconj(P) if x != pari("y")][0]

        def sig(x):
            return pari.subst(pari.lift(x), "y", s) * pari.Mod(1, P)

        t = 2 if str(pari.polredabs(P)) == "y^2 - 2" else 5
        cm = [D for D in QUAD if D < 0 and all(a0[p] == 0 for p in good if kron(D, p) == -1)]
        tw = [D for D in QUAD if all(sig(a0[p]) == kron(D, p) * a0[p] for p in good)]
        certified = False
        for D, kind in [(D, "CM") for D in cm] + [(D, "twist") for D in tw]:
            assert all(M % p == 0 for p in (2, 3, 5) if D % p == 0)
            L = int(pari.lcm(M, D * D))
            B = sturm_bound(L)
            a = pari.mfcoefs(F, B)
            if kind == "CM":
                ok = all(a[n] == kron(D, n) * a[n] for n in range(1, B + 1) if math.gcd(n, abs(D)) == 1)
                note = "CM => elliptic building block"
            else:
                ok = all(sig(a[n]) == kron(D, n) * a[n] for n in range(1, B + 1) if math.gcd(n, abs(D)) == 1)
                hil = [int(pari.hilbert(t, D, v)) for v in (0, 2, 3, 5)]
                assert hil == [1, 1, 1, 1]
                note = f"(t,D)=({t},{D}); exact Hilbert symbols at oo,2,3,5 = {hil}"
            nchk = sum(1 for n in range(1, B + 1) if math.gcd(n, abs(D)) == 1)
            print(f"(V3) M={M:5d} K=Q(sqrt{t}) {kind:5s} D={D:4d} L={L:5d} B(L)={B:5d} "
                  f"coefficients checked={nchk:4d} verified={ok}  {note}")
            certified |= ok
        assert certified, (M, str(P))

    # (V4)
    PR = [p for p in range(7, 400) if pari.isprime(p)]
    x = pari("x")
    triv = []
    for M, chi, mf in spaces():
        if chi != 1:
            continue
        for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
            a = pari.mfcoefs(F, 400)
            d = int(pari.poldegree(P))
            triv.append((d, [pari.charpoly(a[p]) if d > 1 else x - a[p] for p in PR]))
    k = {1: 0, 2: 1, 4: 2, 3: 3}  # discrete log base 2 modulo 5
    n5 = 0
    for M, chi, mf in spaces():
        if chi != 5:
            continue
        for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
            n5 += 1
            d = int(pari.poldegree(P))
            a = pari.mfcoefs(F, 400)
            I = pari.Mod(pari.lift(pari.nfroots(P, pari("x^2+1"))[0]), P)
            cps = [pari.charpoly(I ** (-k[p % 5]) * a[p]) for p in PR]
            assert any(d % dh == 0 and all(cps[i] == ch[i] ** (d // dh) for i in range(len(PR)))
                       for dh, ch in triv), (M, P)
    print(f"(V4) all {n5} chi_5-orbits are quartic twists of trivial-character orbits (7 <= p < 400)")

    # (V5)
    PT = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    tr = {p: 0 for p in PT}
    for M, chi, mf in spaces():
        m = sigma0(N // M)
        for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
            a = pari.mfcoefs(F, 43)
            for p in PT:
                tr[p] += m * (int(pari.trace(a[p])) if int(pari.poldegree(P)) > 1 else int(a[p]))
    print("(V5) traces of Frobenius on H^1:", [tr[p] for p in PT])
    assert [tr[p] for p in PT] == [-8, -84, -18, -30, -12, 8, -34, -64, 6, -246, -36]
    print("All checks passed.")


if __name__ == "__main__":
    main()
