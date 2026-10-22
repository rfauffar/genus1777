# Stage 8: for every block: (i) rank witness: the candidate coordinate matrix (rows = 354 primes, columns = (orbit, y-power))
# reduced mod a prime ell = 1 (mod 24) (zeta_24 -> fixed primitive root of unity) has full column rank -> record a nonsingular
# square minor (row/column indices, determinant mod ell); (ii) the unique multiplicities n_orbit; (iii) the EXACT identity
#   sum_orbit n_orbit * Tr_{K_f/Q(chi)} a_p(f) = -N_p(rho)   in Q(zeta_24) for all 354 primes.
import pickle, os, cypari2
from fractions import Fraction as Fr
import cyclo
pari=cypari2.Pari(); pari.allocatemem(int(2*10**9))
R=pickle.load(open('traces.pkl','rb')); plan=pickle.load(open('plan.pkl','rb')); PR=R['PR']
ell=10**9
while True:
    ell=int(pari(f'nextprime({ell}+1)'))
    if ell%24==1: break
Z=int(pari(f'lift(znprimroot({ell})^(({ell}-1)/24))'))
def modl(v):
    s=0
    for j,c in enumerate(v):
        c=Fr(c); assert c.denominator%ell
        s=(s+c.numerator*pow(c.denominator,-1,ell)*pow(Z,j,ell))%ell
    return s
out={}
for bi,r in enumerate(R['reps']):
    orbs=[]
    for (M,ks) in plan[bi]:
        s=pickle.load(open(f'spaces/{M}_{"_".join(map(str,ks))}.pkl','rb'))
        for o in s['orbits']: orbs.append(((M,tuple(ks)),o))
    cols=[[modl(o['coords'][ip][j]) for ip in range(len(PR))] for (_,o) in orbs for j in range(o['d'])]
    nc=len(cols)
    target=[tuple(-c for c in R['N'][p][bi]) for p in PR]
    if nc==0:
        ok=all(all(c==0 for c in t) for t in target); out[bi]=dict(rep=r,ncols=0,mult={},identity=ok); continue
    def pivot_rows(colsT,k):
        # choose rows (primes) one by one, keeping those independent of the previous ones, until k are found
        basis=[]; piv=[]; rows=[]
        for i in range(len(PR)):
            v=[colsT[j][i] for j in range(len(colsT))]
            for b,pc in zip(basis,piv):
                if v[pc]:
                    f=v[pc]; v=[(x-f*y)%ell for x,y in zip(v,b)]
            pc=next((j for j,x in enumerate(v) if x),None)
            if pc is None: continue
            inv=pow(v[pc],-1,ell); v=[(x*inv)%ell for x in v]
            basis.append(v); piv.append(pc); rows.append(i)
            if len(rows)==k: break
        return rows
    def det_mod(M):
        M=[row[:] for row in M]; n_=len(M); d=1
        for c in range(n_):
            p_=next((i for i in range(c,n_) if M[i][c]%ell),None)
            if p_ is None: return 0
            if p_!=c: M[c],M[p_]=M[p_],M[c]; d=-d
            d=d*M[c][c]%ell; inv=pow(M[c][c],-1,ell)
            for i in range(c+1,n_):
                f=M[i][c]*inv%ell
                if f: M[i]=[(x-f*y)%ell for x,y in zip(M[i],M[c])]
        return d%ell
    def solve_mod(Msq,b):
        n_=len(Msq); Aug=[Msq[i][:]+[b[i]] for i in range(n_)]
        for c in range(n_):
            p_=next(i for i in range(c,n_) if Aug[i][c]%ell)
            Aug[c],Aug[p_]=Aug[p_],Aug[c]; inv=pow(Aug[c][c],-1,ell); Aug[c]=[(x*inv)%ell for x in Aug[c]]
            for i in range(n_):
                if i!=c and Aug[i][c]:
                    f=Aug[i][c]; Aug[i]=[(x-f*y)%ell for x,y in zip(Aug[i],Aug[c])]
        return [Aug[i][n_] for i in range(n_)]
    rows=pivot_rows(cols,nc); rk=len(rows)
    det=det_mod([[cols[j][i] for j in range(nc)] for i in rows]) if rk==nc else None
    cl=list(range(nc))
    tcols=[[modl(o['traces'][ip]) for ip in range(len(PR))] for (_,o) in orbs]
    trows=pivot_rows(tcols,len(orbs))
    bvec=[modl(target[i]) for i in range(len(PR))]
    n=solve_mod([[tcols[j][i] for j in range(len(orbs))] for i in trows],[bvec[i] for i in trows]) if len(trows)==len(orbs) else [None]*len(orbs)
    n=[(x if x<ell//2 else x-ell) if x is not None else -1 for x in n]
    ok=all(x>=0 for x in n)
    if ok:
        for ip,p in enumerate(PR):
            s=[Fr(0)]*cyclo.DEG
            for (_,o),k in zip(orbs,n):
                if k:
                    for j in range(cyclo.DEG): s[j]+=k*o['traces'][ip][j]
            if tuple(s)!=tuple(Fr(c) for c in target[ip]): ok=False; break
    mult={(sp[0],sp[1],o['index']):k for (sp,o),k in zip(orbs,n) if k}
    out[bi]=dict(rep=r,ncols=nc,rank=rk,full_rank=(rk==nc),witness=dict(rows=rows,det_mod_ell=det),mult=mult,identity=ok)
    print(f"block {bi:2d} rep {r}: columns {nc:3d} rank {rk:3d} {'FULL' if rk==nc else 'DEFICIENT'}  exact identity {'OK' if ok else 'FAIL'}  orbits used {len(mult)}",flush=True)
pickle.dump(dict(ell=ell,Z=Z,blocks=out),open('match.pkl','wb'))
print("ell =",ell,"; blocks with full rank and exact identity:",sum(1 for v in out.values() if v.get('identity') and v.get('full_rank',v['ncols']==0)),"of",len(out))
