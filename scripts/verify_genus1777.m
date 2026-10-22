// verify_genus1777.m
//
// Independent Magma check of computations (C1) and (C2) of
//   "A modular curve of genus 1777 with completely decomposable Jacobian".
//
// NOTE: this script was written for referees who have Magma; it has NOT been run by
// the author, who verified (C1) and (C2) with the PARI/GP scripts verify_genus1777.py
// and crosscheck_modsym.py.  Please report any syntax problem.
//
// (C1) For every M | 4800: decompose S_2^new(Gamma_0(M)) into Galois orbits of newforms,
//      check that the orbit degrees add up to the dimension, and that
//      sum_M sigma_0(4800/M) dim S_2^new(M) = 889; collect the non-rational orbits
//      (expected: 16, with fields Q(sqrt 2) and Q(sqrt 5)).
// (C2) For each non-rational orbit, find a self-twist (CM) or inner twist (sigma, chi_D)
//      with D supported on {2,3,5}, and prove it by comparing coefficients through the
//      Sturm bound B(L), L = lcm(M, D^2), as in Lemma 3.5.

N  := 4800;
Ds := [-3, -4, 5, -8, 8, 12, -15, -20, 24, -24, 40, -40, 60, 120, -120];
SturmB := func< L | Ceiling( (L/6) * &*[ 1 + 1/p : p in PrimeDivisors(L) ] ) >;

total  := 0;
nonrat := [];
for M in Divisors(N) do
    S := NewSubspace(CuspForms(Gamma0(M), 2));
    d := Dimension(S);
    if d eq 0 then continue; end if;
    orbs := Newforms(S);
    degs := [ #O : O in orbs ];
    assert &+degs eq d;                      // the list of orbits is complete
    total +:= NumberOfDivisors(N div M) * d;
    for O in orbs do
        if #O gt 1 then Append(~nonrat, <M, O[1]>); end if;
    end for;
    printf "M = %o: dim S2new = %o, orbit degrees %o\n", M, d, Sort(degs);
end for;
printf "weighted total = %o (expected 889); non-rational orbits = %o (expected 16)\n",
       total, #nonrat;
assert total eq 889 and #nonrat eq 16;

for t in nonrat do
    M := t[1];  f := t[2];
    K := BaseRing(Parent(f));
    sig := [ s : s in Automorphisms(K) | s(K.1) ne K.1 ][1];
    good := [ p : p in PrimesUpTo(400) | M mod p ne 0 ];
    a := [ K!Coefficient(f, n) : n in [1..400] ];
    certified := false;
    for D in Ds do
        if not &and[ M mod p eq 0 : p in PrimeDivisors(Abs(D)) ] then continue; end if;
        istwist := &and[ sig(a[p]) eq KroneckerSymbol(D, p) * a[p] : p in good ];
        iscm    := D lt 0 and &and[ a[p] eq 0 : p in good | KroneckerSymbol(D, p) eq -1 ];
        if not (istwist or iscm) then continue; end if;
        L := LCM(M, D^2);  B := SturmB(L);
        c := [ K!Coefficient(f, n) : n in [1..B] ];
        if istwist then
            ok := &and[ sig(c[n]) eq KroneckerSymbol(D, n) * c[n] : n in [1..B] | GCD(n, D) eq 1 ];
            printf "M = %o, field %o: inner twist chi_%o, L = %o, B(L) = %o, verified: %o\n",
                   M, DefiningPolynomial(K), D, L, B, ok;
            certified := certified or ok;
        end if;
        if iscm then
            ok := &and[ c[n] eq KroneckerSymbol(D, n) * c[n] : n in [1..B] | GCD(n, D) eq 1 ];
            printf "M = %o, field %o: CM by D = %o, L = %o, B(L) = %o, verified: %o\n",
                   M, DefiningPolynomial(K), D, L, B, ok;
            certified := certified or ok;
        end if;
    end for;
    assert certified;
end for;
print "All Magma checks passed.";
