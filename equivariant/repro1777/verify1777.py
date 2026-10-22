#!/usr/bin/env python3
# Standalone verifier for the genus-1777 certificate. Pure Python 3 (standard library only): no GAP, no PARI.
# Checks: (V1) genus by coset enumeration; (V2) block dimensions sum to the genus; (V3) the exact identity
#   sum_f n_f Tr a_p(f) = -N_p(rho) in Q(zeta_24) for every block and all 354 primes; (V4) rank witnesses: the candidate
#   matrix restricted to the witness rows has nonzero determinant mod ell (zeta_24 -> given root of unity); (V5) ellipticity
#   of every orbit used: CM identities / inner-twist identities up to the recorded Sturm bounds, field automorphisms,
#   and Quer's Brauer class via recomputed square classes and Hilbert symbols.
import gzip, json, math, sys, re
from fractions import Fraction as Fr
cert=json.load(gzip.open(sys.argv[1] if len(sys.argv)>1 else 'cert1777_full.json.gz','rt'))
ok_all=True
def report(name,ok,extra=''):
    global ok_all; ok_all&=ok; print(f"[{'PASS' if ok else 'FAIL'}] {name} {extra}")
# ---------- V1: coset enumeration ----------
N=360
def m4(x,y): a,b,c,d=x; e,f,g,h=y; return ((a*e+b*g)%N,(a*f+b*h)%N,(c*e+d*g)%N,(c*f+d*h)%N)
def crt(x8,x9,x5): return tuple((v8*45*pow(45,-1,8)+v9*40*pow(40,-1,9)+v5*72*pow(72,-1,5))%N for v8,v9,v5 in zip(x8,x9,x5))
g=cert['group']; Gam=set(crt(a,b,c) for a in map(tuple,g['G8']) for b in map(tuple,g['G9']) for c in map(tuple,g['G5']))
ok=all(m4(x,y) in Gam for x in list(Gam)[:200] for y in list(Gam)[:200])
def canon(x): return min(m4(h,x) for h in Gam)
S=(0,N-1,1,0); T=(1,1,0,1); ST=m4(S,T)
st=canon((1,0,0,1)); idx={st:0}; reps=[st]; q=0
while q<len(reps):
    for Mx in (S,T):
        y=canon(m4(reps[q],Mx))
        if y not in idx: idx[y]=len(reps); reps.append(y)
    q+=1
n=len(reps); pT=[idx[canon(m4(x,T))] for x in reps]; pS=[idx[canon(m4(x,S))] for x in reps]; pST=[idx[canon(m4(x,ST))] for x in reps]
e2=sum(pS[i]==i for i in range(n)); e3=sum(pST[i]==i for i in range(n)); seen=[0]*n; c=0
for i in range(n):
    if not seen[i]:
        c+=1; j=i
        while not seen[j]: seen[j]=1; j=pT[j]
genus=1+Fr(n,12)-Fr(e2,4)-Fr(e3,3)-Fr(c,2)
report("V1 genus by coset enumeration",ok and genus==1777,f"(index {n}, e2 {e2}, e3 {e3}, cusps {c}, genus {genus})")
# ---------- V2 ----------
tot=sum(Fr(b['dim']) for b in cert['blocks'])
report("V2 block dimensions sum to the genus",tot==1777,f"({len(cert['blocks'])} blocks, total {tot})")
# ---------- V3: exact identities ----------
def V(v): return [Fr(x) for x in v]
bad=[]
for b in cert['blocks']:
    for ip in range(len(cert['primes'])):
        s=[Fr(0)]*8
        for M,ks,oi,mult in b['multiplicities']:
            sp=cert['spaces'][f"{M}_{'_'.join(map(str,ks))}"]; o=[x for x in sp['orbits'] if x['index']==oi][0]
            s=[x+mult*y for x,y in zip(s,V(o['traces'][ip]))]
        if s!=V(b['minus_N_p'][ip]): bad.append((b['index'],cert['primes'][ip])); break
report("V3 exact identity sum n_f Tr a_p(f) = -N_p(rho) in Q(zeta_24), all blocks x all primes",not bad,str(bad[:3]))
# ---------- V4: rank witnesses ----------
ell=cert['ell']; Z=cert['zeta24_mod_ell']
def modl(v):
    s=0
    for j,x in enumerate(v):
        x=Fr(x); s=(s+x.numerator*pow(x.denominator,-1,ell)*pow(Z,j,ell))%ell
    return s
def det_mod(M):
    M=[r[:] for r in M]; n_=len(M); d=1
    for col in range(n_):
        p_=next((i for i in range(col,n_) if M[i][col]),None)
        if p_ is None: return 0
        if p_!=col: M[col],M[p_]=M[p_],M[col]; d=-d
        d=d*M[col][col]%ell; inv=pow(M[col][col],-1,ell)
        for i in range(col+1,n_):
            f=M[i][col]*inv%ell
            if f: M[i]=[(x-f*y)%ell for x,y in zip(M[i],M[col])]
    return d%ell
assert pow(Z,24,ell)==1 and pow(Z,12,ell)!=1 and pow(Z,8,ell)!=1
fails=[]
for b in cert['blocks']:
    cols=[]
    for M,ks in b['candidate_spaces']:
        sp=cert['spaces'][f"{M}_{'_'.join(map(str,ks))}"]
        for o in sp['orbits']:
            for j in range(o['rel_degree']): cols.append([modl(o['coords_at_rows'][str(i)][j]) for i in b['rank_witness_rows']])
    if not cols: continue
    nr=len(b['rank_witness_rows'])
    if nr!=len(cols): fails.append((b['index'],'size')); continue
    d=det_mod([[cols[j][i] for j in range(len(cols))] for i in range(nr)])
    if d==0 or d!=b['det_mod_ell']: fails.append((b['index'],d))
report("V4 rank witnesses (nonzero minors mod ell; candidate columns linearly independent)",not fails,str(fails[:3]))
# ---------- V5: ellipticity ----------
def parse_poly(s):
    # polynomial in y with rational coefficients -> dict deg->Fraction
    s=s.replace(' ','')
    if s=='' : return {}
    terms=re.findall(r'[+-]?[^+-]+',s); P={}
    for t in terms:
        coef=Fr(1); deg=0
        sign=-1 if t.startswith('-') else 1; t=t.lstrip('+-')
        if 'y' in t:
            pre,post=t.split('y',1)
            pre=pre.rstrip('*'); coef=Fr(pre) if pre not in ('',) else Fr(1)
            deg=int(post[1:]) if post.startswith('^') else 1
        else: coef=Fr(t)
        P[deg]=P.get(deg,0)+sign*coef
    return {k:v for k,v in P.items() if v}
def pmod(P,F):
    P=dict(P); d=max(F)
    while P and max(P)>=d:
        k=max(P); c=P.pop(k)/F[d]
        for j,v in F.items():
            if j!=d: P[k-d+j]=P.get(k-d+j,0)-c*v
        P={a:b for a,b in P.items() if b}
    return P
def pmul(A,B):
    C={}
    for i,x in A.items():
        for j,y in B.items(): C[i+j]=C.get(i+j,0)+x*y
    return {k:v for k,v in C.items() if v}
def padd(A,B,s=1):
    C=dict(A)
    for k,v in B.items(): C[k]=C.get(k,0)+s*v
    return {k:v for k,v in C.items() if v}
def subst(P,Sg,F):
    R={}; pw={0:Fr(1)}
    for k in range(max(P)+1 if P else 0):
        if k in P: R=padd(R,{a:b*P[k] for a,b in pw.items()})
        pw=pmod(pmul(pw,Sg),F)
    return pmod(R,F)
def kron(D,n):
    # Kronecker symbol (D/n) for n>0
    if math.gcd(D,n)!=1: return 0
    r=1
    while n%2==0:
        n//=2
        if D%8 in (3,5): r=-r
    a=D%n; m=n; t=1
    while a:
        while a%2==0:
            a//=2
            if m%8 in (3,5): t=-t
        a,m=m,a
        if a%4==3 and m%4==3: t=-t
        a%=m
    return r*t if m==1 else 0
def sturm(L):
    pr=[p for p in range(2,L+1) if L%p==0 and all(p%q for q in range(2,int(p**.5)+1))]
    return math.ceil(Fr(L)*math.prod(Fr(p+1,p) for p in pr)/6)
def core(x):
    x=int(x); s=-1 if x<0 else 1; x=abs(x); r=1; p=2
    while p*p<=x:
        while x%(p*p)==0: x//=p*p
        if x%p==0: r*=p; x//=p
        p+=1
    return s*r*x
def hilbert(a,b,p):
    # Hilbert symbol (a,b)_p for nonzero integers; p=0 means the real place
    if p==0: return -1 if (a<0 and b<0) else 1
    def split(x):
        v=0
        while x%p==0: x//=p; v+=1
        return v,x
    va,ua=split(a); vb,ub=split(b)
    if p!=2:
        leg=lambda u: pow(u%p,(p-1)//2,p)
        e=(-1)**(va*vb*((p-1)//2)%2) if False else 1
        s=pow(-1,va*vb*((p-1)//2))*(1 if leg(ua)**vb%p==1 or vb%2==0 else -1)*(1 if leg(ub)**va%p==1 or va%2==0 else -1)
        lu=lambda u: 1 if leg(u)==1 else -1
        return pow(-1,va*vb*((p-1)//2))*(lu(ua)**vb)*(lu(ub)**va)
    eps=lambda u: ((u-1)//2)%2; omg=lambda u: ((u*u-1)//8)%2
    return (-1)**((eps(ua)*eps(ub)+va*omg(ub)+vb*omg(ua))%2)
nfail=[]
for k,v in cert['orbits'].items():
    if v['type']=='rational': continue
    F=parse_poly(v['field_abs']) if 'y' in v['field_abs'] else {1:Fr(1),0:Fr(0)}
    a=[pmod(parse_poly(x) if 'y' in x else ({0:Fr(x)} if Fr(x) else {}),F) for x in v['a_n']]
    M=v['M']
    if v['type']=='CM':
        D=v['D']; B=v['B']
        good=B>=sturm(math.lcm(M,D*D)) and all(not a[n] or kron(D,n)==1 for n in range(1,B+1) if math.gcd(n,M*abs(D))==1)
        if not good: nfail.append((k,'CM'))
        continue
    deg=max(F); autos=[t['sigma'] for t in v['twists']]
    good=len(autos)==deg
    for t in v['twists']:
        sg=parse_poly(t['sigma'])
        if subst(F,sg,F): good=False                                  # sigma(y) must be a root of F
        if t['sigma']=='y': continue
        D=t['D']; B=t['checked_up_to']
        if B<(sturm(math.lcm(M,D*D)) if D!=1 else sturm(M)): good=False
        for nn in range(1,B+1):
            if math.gcd(nn,M*abs(D))!=1: continue
            if padd(subst(a[nn],sg,F),{e:c*kron(D,nn) for e,c in a[nn].items()},-1): good=False; break
    # Quer: recompute square classes of b_p^2 = a_p^2/eps(p) for 7<=p<=2400, p !| M, and the Hilbert symbols
    if good and v.get('hilbert_symbols') is not None:
        # recompute the square classes of b_p^2 = a_p^2/eps(p) (eps(p) = chi_5-type sign recorded via kronecker data of the twists)
        # nebentypus eps from the stored character of the orbit's space (generators, orders, exponents)
        spc=cert['spaces'][v['space']]; Mi=spc['M']; tabe={}
        import itertools as _it
        for e_ in _it.product(*[range(c_) for c_ in spc['cyc']]):
            u_=1
            for g_,x_ in zip(spc['gens'],e_): u_=u_*pow(g_,x_,Mi)%Mi
            tabe[u_]=sum(Fr(k_*x_,c_) for k_,x_,c_ in zip(spc['ks'],e_,spc['cyc']))%1
        assert all(val in (0,Fr(1,2)) for val in tabe.values()), "nebentypus must be real here"
        epsf=lambda p_: -1 if tabe[p_%Mi]==Fr(1,2) else 1
        sq={}
        for p in [x for x in range(7,min(len(a)-1,2400)+1) if all(x%r for r in range(2,int(x**.5)+1))]:
            if M%p==0 or not a[p]: continue
            b2=pmod(pmul(a[p],a[p]),F)
            if any(k_!=0 for k_ in b2): good=False; break                 # a_p^2 must be rational
            val=b2.get(0,0)*epsf(p)
            sq[p]=core(val)
        gens=[]; span={1:0}
        for x in sorted(set(sq.values()),key=abs):
            if x not in span:
                kk=len(gens); gens.append(x); span.update({core(s_*x):(m_|(1<<kk)) for s_,m_ in list(span.items())})
        if sorted(gens)!=sorted(t_ for t_,D_ in v['quer_pairs']): good=False
        for kk,(t_,D_) in enumerate(v['quer_pairs']):
            if t_ not in gens: good=False; break
            kpos=gens.index(t_)
            if not all(((span[sq[p]]>>kpos)&1)==(1 if kron(D_,p)==-1 else 0) for p in sq): good=False
        sym_ok=all(math.prod(hilbert(t_,D_,int(q)) for t_,D_ in v['quer_pairs'])==1 for q in v['hilbert_symbols'])
        if not sym_ok or not v['brauer_split']: good=False
    if not good: nfail.append((k,'twist/Brauer'))
report("V5 ellipticity of every non-rational orbit (CM / inner twists / split Brauer class)",not nfail,str(nfail[:3]))
print("\nALL CHECKS PASSED" if ok_all else "\nSOME CHECKS FAILED")
