# Stage 3: the group Gamma-bar in SL2(Z/360) = G8 x G9 x G5, its genus by coset enumeration, local invariant dimensions,
# Chevalley-Weil multiplicities g_rho (exact), and the representation blocks of Jac(X_Gamma).
import itertools, math, pickle
from fractions import Fraction
import cyclo
C=pickle.load(open('classes.pkl','rb')); T=pickle.load(open('tables.pkl','rb'))
N=360; MODS=(8,9,5)
def m4(x,y,n): a,b,c,d=x; e,f,g,h=y; return ((a*e+b*g)%n,(a*f+b*h)%n,(c*e+d*g)%n,(c*f+d*h)%n)
# --- local factors (upper-unipotent version mod 3; conjugate to the lower-unipotent one) ---
G8=[(l,0,0,l) for l in range(8) if (l*l)%8==1]                                     # scalars lambda, lambda^2=1
G9=[x for x in itertools.product(range(9),repeat=4) if (x[0]*x[3]-x[1]*x[2])%9==1
    and x[2]%3==0 and ((x[0]%3==1 and x[3]%3==1) or (x[0]%3==2 and x[3]%3==2))]     # +-(1 b;0 1) mod 3
G5=[(1,0,0,1),(4,0,0,4)]                                                           # +-I mod 5
print("local orders:",len(G8),len(G9),len(G5))
# --- genus by coset enumeration of PSL2: cosets of Gamma-bar in SL2(Z/360) ---
def crt(x8,x9,x5):
    return tuple((v8*45*pow(45,-1,8)+v9*40*pow(40,-1,9)+v5*72*pow(72,-1,5))%N for v8,v9,v5 in zip(x8,x9,x5))
Gam=set(crt(a,b,c) for a in G8 for b in G9 for c in G5)
def canon(g): return min(m4(h,g,N) for h in Gam)
S=(0,N-1,1,0); Tm=(1,1,0,1); ST=m4(S,Tm,N)
start=canon((1,0,0,1)); idx={start:0}; reps=[start]; q=0
while q<len(reps):
    for M in (S,Tm):
        y=canon(m4(reps[q],M,N))
        if y not in idx: idx[y]=len(reps); reps.append(y)
    q+=1
n=len(reps)
pT=[idx[canon(m4(x,Tm,N))] for x in reps]; pS=[idx[canon(m4(x,S,N))] for x in reps]; pST=[idx[canon(m4(x,ST,N))] for x in reps]
e2=sum(pS[i]==i for i in range(n)); e3=sum(pST[i]==i for i in range(n))
seen=[False]*n; cusps=0
for i in range(n):
    if not seen[i]:
        cusps+=1; j=i
        while not seen[j]: seen[j]=True; j=pT[j]
genus=1+Fraction(n,12)-Fraction(e2,4)-Fraction(e3,3)-Fraction(cusps,2)
print("coset enumeration: index",n,"e2",e2,"e3",e3,"cusps",cusps,"genus",genus)
# --- exact characters ---
def chi(m,i,x):
    cid,sizes,order=C[m]; return cyclo.from_coeffs(T[m]['ex'][i][cid[x]],T[m]['L'])
def inv_dim_local(m,i,H):
    s=cyclo.ZERO
    for h in H: s=cyclo.add(s,chi(m,i,h))
    assert cyclo.is_rational(s) and s[0]%len(H)==0
    return s[0]//len(H)
dloc={m:[inv_dim_local(m,i,G) for i in range(len(T[m]['ex']))] for m,G in ((8,G8),(9,G9),(5,G5))}
# --- Chevalley-Weil at level 360: g_rho for rho=(i8,i9,i5) with dim rho^Gamma>0 ---
def gen(gens):
    Sx={(1,0,0,1)}; fr=[(1,0,0,1)]
    while fr:
        nf=[]
        for x in fr:
            for g in gens:
                y=m4(x,g,N)
                if y not in Sx: Sx.add(y); nf.append(y)
        fr=nf
    return sorted(Sx)
mI=(N-1,0,0,N-1)
H={'pm1':gen([mI]),'C4':gen([S]),'C6':gen([ST,mI]),'U':gen([Tm,mI])}
def red(x,m): return tuple(v%m for v in x)
def inv_dim(rho,Hs):
    s=cyclo.ZERO
    for h in Hs:
        s=cyclo.add(s,cyclo.mul(cyclo.mul(chi(8,rho[0],red(h,8)),chi(9,rho[1],red(h,9))),chi(5,rho[2],red(h,5))))
    assert cyclo.is_rational(s) and s[0]%len(Hs)==0
    return s[0]//len(Hs)
def dim(m,i):
    v=cyclo.from_coeffs(T[m]['ex'][i][C[m][0][(1,0,0,1)]],T[m]['L']); assert cyclo.is_rational(v); return v[0]   # GAP coefficients are not unique; reduce first
def is_det_char(m,i):
    cid=C[m][0]
    return dim(m,i)==1 and all(cyclo.from_coeffs(T[m]['ex'][i][cid[x]],T[m]['L'])==cyclo.ONE for x in cid if (x[0]*x[3]-x[1]*x[2])%m==1)
support=[(a,b,c) for a in range(len(dloc[8])) for b in range(len(dloc[9])) for c in range(len(dloc[5])) if dloc[8][a]*dloc[9][b]*dloc[5][c]>0]
print("representations with Gamma-invariants:",len(support))
g={}
for rho in support:
    det=is_det_char(8,rho[0]) and is_det_char(9,rho[1]) and is_det_char(5,rho[2])
    val=inv_dim(rho,H['pm1'])-inv_dim(rho,H['C4'])-inv_dim(rho,H['C6'])-inv_dim(rho,H['U'])+(2 if det else 0)
    assert val%2==0; g[rho]=val//2
phiN=sum(1 for u in range(N) if math.gcd(u,N)==1)
contrib={rho:dloc[8][rho[0]]*dloc[9][rho[1]]*dloc[5][rho[2]]*g[rho] for rho in support if g[rho]>0}
tot=Fraction(sum(contrib.values()),phiN)
print("genus from representation theory: sum dim(rho^Gamma) g_rho / phi(360) =",tot)
pickle.dump(dict(G8=G8,G9=G9,G5=G5,coset=dict(index=n,e2=e2,e3=e3,cusps=cusps,genus=str(genus)),dloc=dloc,g=g,contrib=contrib,phiN=phiN),open('group.pkl','wb'))
