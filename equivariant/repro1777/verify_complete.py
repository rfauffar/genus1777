#!/usr/bin/env python3
# V7: completeness of the candidate newform sets, recomputed from the GAP tables (pure Python):
#  (a) local conductor bounds of each block representative via the corrected lemma (H_{k,n}-invariants, n>=1; spherical n=0);
#  (b) its central character omega on (Z/360)^*;
#  (c) the stored candidate spaces are exactly {(M, eps): M = 2^a 3^b 5^c within the bounds, M >= 11, eps in {omega, conj(omega)}
#      factoring through (Z/M)^*}, with each stored character (generators, orders, exponents) evaluated and compared;
#  (d) each newspace dimension equals the Cohen-Oesterle formula (via the newform decomposition with beta(p)=-2, beta(p^2)=1),
#      and equals the number of newforms listed (sum of relative degrees over Q(chi) of the orbits).
import gzip, json, sys, math, itertools, cmath
from fractions import Fraction as Fr
from vcommon import *
cert=json.load(gzip.open(sys.argv[1] if len(sys.argv)>1 else 'cert1777_full.json.gz','rt'))
TB=load_tables(cert['gap_tables']); N=360; ok_all=True
def inv_dim(m,i,Hs):
    cnt={}
    for x in Hs: k=TB[m]['cid'][x]; cnt[k]=cnt.get(k,0)+1
    s=[0]*8
    for k,c in cnt.items():
        v=TB[m]['vals'][i][k]
        for j in range(8): s[j]+=c*v[j]
    assert all(x==0 for x in s[1:]) and s[0]%len(Hs)==0
    return s[0]//len(Hs)
def Hset(m,p,e,k,n):
    S=set()
    for a,b,c,d in itertools.product(range(m),repeat=4):
        if n>0 and (d-1)%(p**min(n,e)): continue
        x=(a,(p**k*b)%m,(p**(n-k)*c)%m if n-k<=e else 0,d)
        if (x[0]*x[3]-x[1]*x[2])%p: S.add(x)
    return S
HS={}
def bound(m,p,e,i):
    for n in range(0,2*e+1):
        for k in range(0,min(n,e)+1):
            if n-k>e: continue
            if (m,k,n) not in HS: HS[(m,k,n)]=Hset(m,p,e,k,n)
            if inv_dim(m,i,HS[(m,k,n)])>0: return n
    return 2*e
def omega_exp(r,u):
    v=ONE; dims=1
    for t,m in enumerate((8,9,5)):
        I=TB[m]['cid'][(1,0,0,1)]; dims*=TB[m]['vals'][r[t]][I][0]
        v=mul24(v,TB[m]['vals'][r[t]][TB[m]['cid'][(u%m,0,0,u%m)]])
    for k in range(24):
        if red24([0]*k+[dims])==v: return k
    raise ValueError
U360=[u for u in range(1,N) if math.gcd(u,N)==1]
def char_table(sp):
    # values (as fractions mod 1) of the stored character on all units mod M, from generators/orders/exponents
    M=sp['M']; tab={}
    for e in itertools.product(*[range(c) for c in sp['cyc']]):
        u=1
        for g,x in zip(sp['gens'],e): u=u*pow(g,x,M)%M
        tab[u]=sum(Fr(k*x,c) for k,x,c in zip(sp['ks'],e,sp['cyc']))%1
    phiM=sum(1 for u in range(1,M) if math.gcd(u,M)==1) if M>1 else 1
    assert len(tab)==phiM, "stored generators do not generate (Z/M)^* freely"
    return tab
def cohen_oesterle(d,epsv):
    # dim S_2(Gamma_0(d), eps) ; epsv: dict unit mod d -> value in Q/Z
    cond=min(c for c in range(1,d+1) if d%c==0 and all(epsv[u]==0 for u in epsv if (u-1)%c==0))
    fac={}; x=d; p=2
    while x>1:
        while x%p==0: fac[p]=fac.get(p,0)+1; x//=p
        p+=1
    mu=d*math.prod(Fr(p+1,p) for p in fac)
    lam=1
    for p,r in fac.items():
        s=0; c=cond
        while c%p==0: s+=1; c//=p
        if 2*s<=r: lam*=(p**(r//2)+p**(r//2-1)) if r%2==0 else 2*p**(r//2)
        else: lam*=2*p**(r-s)
    ev=lambda u: cmath.exp(2j*math.pi*float(epsv[u]))
    s4=sum(ev(x) for x in range(d) if (x*x+1)%d==0 and math.gcd(x,d)==1)
    s3=sum(ev(x) for x in range(d) if (x*x+x+1)%d==0 and math.gcd(x,d)==1)
    val=complex(mu/12)-lam/2-s4/4-s3/3+(1 if all(v==0 for v in epsv.values()) else 0)
    assert abs(val.imag)<1e-9 and abs(val.real-round(val.real))<1e-9
    return round(val.real),cond
def newdim(M,eps360):
    # eps360: dict unit mod 360 -> Q/Z ; returns dim S_2^new(M, eps) (eps must factor mod M)
    def restrict(d):
        tab={}
        for u in U360:
            if u%d in tab and tab[u%d]!=eps360[u]: return None
            tab[u%d]=eps360[u]
        # units mod d not hit by units mod 360 (only if d has primes outside {2,3,5}: not the case here)
        return tab
    beta=lambda n: math.prod({1:-2,2:1}.get(v,0) for v in [sum(1 for _ in iter(lambda: None,None)) for _ in []]) if False else None
    def beta(n):
        r=1; p=2
        while n>1:
            v=0
            while n%p==0: n//=p; v+=1
            r*={0:1,1:-2,2:1}.get(v,0); p+=1
        return r
    tot=0
    for d in [d for d in range(1,M+1) if M%d==0]:
        tab=restrict(d)
        if tab is None: continue
        if d<3:   # trivial dims at tiny levels
            continue
        dim_d,_=cohen_oesterle(d,tab); tot+=beta(M//d)*dim_d
    return tot
fails=[]
for b in cert['blocks']:
    r=tuple(b['rep'])
    bd=(bound(8,2,3,r[0]),bound(9,3,2,r[1]),bound(5,5,1,r[2]))
    om={u:Fr(omega_exp(r,u),24) for u in U360}
    expected=[]
    for a in range(bd[0]+1):
        for bb in range(bd[1]+1):
            for c in range(bd[2]+1):
                M=2**a*3**bb*5**c
                if M<11: continue
                for sgn in (1,-1):
                    eps={u:(sgn*om[u])%1 for u in U360}
                    if all(eps[u]==0 for u in U360 if (u-1)%M==0): expected.append((M,tuple(sorted(eps.items()))))
    expected=set(expected)
    got=set()
    for M,ks in b['candidate_spaces']:
        sp=cert['spaces'][f"{M}_{'_'.join(map(str,ks))}"]; tab=char_table(sp)
        got.add((M,tuple(sorted((u,tab[u%M]) for u in U360))))
    if expected!=got: fails.append((b['index'],'space list',len(expected),len(got))); continue
fails2=[]
for key,sp in cert['spaces'].items():
    tab=char_table(sp)
    eps360={u:tab[u%sp['M']] for u in U360}
    nd=newdim(sp['M'],eps360)
    ncount=sum(o['rel_degree'] for o in sp['orbits'])
    if not (nd==sp['newspace_dim']==ncount): fails2.append((key,nd,sp['newspace_dim'],ncount))
print(f"[{'PASS' if not fails else 'FAIL'}] V7a candidate spaces = all spaces allowed by the recomputed conductor bounds and central characters {fails[:3]}")
print(f"[{'PASS' if not fails2 else 'FAIL'}] V7b newspace dimensions: Cohen-Oesterle = recorded dimension = number of newforms listed ({len(cert['spaces'])} spaces) {fails2[:3]}")
