#!/usr/bin/env python3
"""
Cross-check of computation (C1) by a second algorithm.

For every M | 4800, the new subspace of the space of weight-2 modular symbols for
Gamma_0(M) (sign +1) is split into Hecke-irreducible pieces with PARI's modular-symbols
package (msinit / msnew / mssplit / mshecke), which is independent of the mf package
used in verify_genus1777.py.  We check:

  (X1) the multiset of piece dimensions equals the multiset of coefficient-field degrees
       of the newform orbits returned by mfinit/mfsplit;
  (X2) for p = 7, 11, 13 the characteristic polynomial of T_p on the whole new subspace
       equals the product, over the newform orbits, of the characteristic polynomials
       of a_p(f) over Q;
  (X3) for each piece of dimension > 1, the characteristic polynomials of T_7, T_11,
       T_13, T_17 on the piece coincide with those of a_7, a_11, a_13, a_17 of a
       non-rational mf orbit of the same level (giving a bijection with Table 1).

Runtime: about 15-20 minutes, dominated by level 4800.
"""
import sys
import time
import cypari2

pari = cypari2.Pari()
pari.allocatemem(8 * 10**9)

N = 4800
DIVS = [d for d in range(1, N + 1) if N % d == 0]
PS = [7, 11, 13]
PS3 = [7, 11, 13, 17]
x = pari("x")

ok_all = True
for M in DIVS:
    t0 = time.time()
    mf = pari(f"mfinit([{M},2,1],0)")
    d = int(pari.mfdim(mf))
    if d == 0:
        continue
    # mf data
    orbs = []
    for F, P in zip(pari.mfeigenbasis(mf), pari.mffields(mf)):
        a = pari.mfcoefs(F, 17)
        deg = int(pari.poldegree(P))
        cps = {p: (pari.charpoly(a[p]) if deg > 1 else x - a[p]) for p in PS3}
        orbs.append((deg, cps))
    # modular symbols data
    pari(f"MS=msinit({M},2,1); NW=msnew(MS); SP=mssplit(MS,NW);")
    dims = sorted(int(v) for v in pari("apply(w->matsize(w[1])[2], SP)"))
    x1 = dims == sorted(o[0] for o in orbs)
    x2 = True
    for p in PS:
        cp_ms = pari(f"charpoly(mshecke(MS,{p},NW))")
        cp_mf = pari(1)
        for o in orbs:
            cp_mf *= o[1][p]
        x2 &= (cp_ms == cp_mf)
    x3 = True
    nonrat = [o for o in orbs if o[0] > 1]
    used = [False] * len(nonrat)
    npieces = int(pari("#SP"))
    for i in range(1, npieces + 1):
        if int(pari(f"matsize(SP[{i}][1])[2]")) == 1:
            continue
        cps = {p: pari(f"charpoly(mshecke(MS,{p},SP[{i}]))") for p in PS3}
        hit = [j for j, o in enumerate(nonrat) if not used[j] and all(o[1][p] == cps[p] for p in PS3)]
        if hit:
            used[hit[0]] = True
        else:
            x3 = False
    x3 &= all(used)
    ok = x1 and x2 and x3
    ok_all &= ok
    print(f"M={M:5d} new dim={d:3d} orbits={len(orbs):3d} non-rational={len(nonrat)} "
          f"X1={x1} X2={x2} X3={x3} ({time.time() - t0:.0f}s)", flush=True)

print("All modular-symbol cross-checks passed." if ok_all else "SOME CHECK FAILED")
sys.exit(0 if ok_all else 1)
