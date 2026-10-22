# Exact arithmetic in Z[zeta_L] (L = 24) with integer coefficient vectors reduced modulo the cyclotomic polynomial.
from functools import lru_cache
L=24
def _phi(n):
    # cyclotomic polynomial Phi_n as integer coefficient list (low degree first)
    import sympy
    x=sympy.symbols('x'); P=sympy.Poly(sympy.cyclotomic_poly(n,x),x)
    return [int(c) for c in reversed(P.all_coeffs())]
PHI=_phi(L); DEG=len(PHI)-1
def reduce(v):
    v=list(v)+[0]*max(0,2*L-len(v))
    for k in range(len(v)-1,DEG-1,-1):
        c=v[k]
        if c:
            for j in range(DEG+1): v[k-DEG+j]-=c*PHI[j]
    return tuple(v[:DEG])
def from_coeffs(coeffs,Lc):
    # element sum c_k zeta_Lc^k  (Lc | L)
    v=[0]*L
    for k,c in enumerate(coeffs): v[(k*(L//Lc))%L]+=c
    return reduce(v)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def mul(a,b):
    v=[0]*(2*DEG)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y: v[i+j]+=x*y
    return reduce(v)
ZERO=tuple([0]*DEG); ONE=tuple([1]+[0]*(DEG-1))
def is_rational(a): return all(x==0 for x in a[1:])
