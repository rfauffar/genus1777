# Stage 7: candidate newform spaces for every block representative rho (complete by the proved conductor bounds of Stage 5
# and the central character of rho), with exact Hecke data from PARI:
#   for each Galois orbit over Q(chi): relative degree d, coordinates of a_p in the basis 1,y,..,y^(d-1) over Q(chi)  (exact),
#   and the relative trace Tr_{K_f/Q(chi)} a_p (exact), for the 354 primes of Stage 6.  Q(chi) is embedded in Q(zeta_24).
import pickle, os, time, math, cmath
from fractions import Fraction as Fr
import cypari2, cyclo
pari=cypari2.Pari(); pari.allocatemem(int(3*10**9)); pari('default(nbthreads,1)')
C=pickle.load(open('classes.pkl','rb')); T=pickle.load(open('tables.pkl','rb')); BD=pickle.load(open('bounds.pkl','rb')); R=pickle.load(open('traces.pkl','rb'))
PR=R['PR']; reps=R['reps']; N=360
chi={m:[[cyclo.from_coeffs(T[m]['ex'][i][k],T[m]['L']) for k in range(len(C[m][1]))] for i in range(len(T[m]['ex']))] for m in (8,9,5)}
Z24=[cmath.exp(2j*math.pi*k/24) for k in range(24)]
def omega_exp(r,u):
    # central character of rho at the scalar u (u coprime to 360): returns k with omega(u) = zeta_24^k (exact)
    v=cyclo.ONE
    for t,m in enumerate((8,9,5)):
        cid=C[m][0]; dim=chi[m][r[t]][cid[(1,0,0,1)]][0]
        w=chi[m][r[t]][cid[(u%m,0,0,u%m)]]
        v=cyclo.mul(v,w); 
    dims=1
    for t,m in enumerate((8,9,5)): dims*=chi[m][r[t]][C[m][0][(1,0,0,1)]][0]
    for k in range(24):
        if cyclo.from_coeffs([0]*k+[dims],24)==v: return k
    raise ValueError("central character value is not a root of unity times dim")
def charvec(M,r,conj):
    # Conrey-free description: exponents ks on znstar(M,1).gen with chi(gen_i)=exp(2 pi i ks_i/cyc_i); None if omega does not factor mod M
    gens=[int(x) for x in pari(f'lift(znstar({M},1).gen)')]; cyc=[int(x) for x in pari(f'znstar({M},1).cyc')]
    ks=[]
    for g,cy in zip(gens,cyc):
        gg=g
        while math.gcd(gg,N)!=1: gg+=M
        k=omega_exp(r,gg); k=(-k)%24 if conj else k
        if (k*cy)%24: return None
        ks.append((k*cy//24)%cy)
    for tt in range(1,60):
        u=1+tt*M
        if math.gcd(u,N)==1 and omega_exp(r,u)%24!=0: return None
    return ks
def to24(poly_t,o):
    # element of Q(zeta_o) given as polynomial in t (t=zeta_o) -> reduced Z[zeta_24]-coefficient vector (Fractions)
    v=[Fr(0)]*24
    for j,c in enumerate(poly_t): v[(j*(24//o))%24]+=c
    w=[Fr(0)]*cyclo.DEG
    # reduce with exact rational arithmetic
    vv=list(v)+[Fr(0)]*24
    for k in range(len(vv)-1,cyclo.DEG-1,-1):
        c=vv[k]
        if c:
            for j in range(cyclo.DEG+1): vv[k-cyclo.DEG+j]-=c*cyclo.PHI[j]
    return tuple(vv[:cyclo.DEG])
def tcoeffs(expr,o):
    a=pari(f'liftall({expr})')
    if o<=1: return [Fr(str(a))]
    return [Fr(str(pari(f'polcoef({a},{j},t)'))) for j in range(o)]
os.makedirs('spaces',exist_ok=True)
def space(M,ks):
    fn=f'spaces/{M}_{"_".join(map(str,ks))}.pkl'
    if os.path.exists(fn): return pickle.load(open(fn,'rb'))
    pari(f'G=znstar({M},1); chi={list(ks)}; mf=mfinit([{M},2,[G,chi]],0); o=charorder(G,chi)')
    o=int(pari('o')); res={'M':M,'ks':list(ks),'o':o,'dim':int(pari('mfdim(mf)')),'orbits':[]}
    if res['dim']>0:
        pari(f'S=mfsplit(mf); vF=S[1]; vK=S[2]; Cm=mfcoefs(mf,{PR[-1]})')
        for i in range(1,int(pari('#vK'))+1):
            pari(f'w=Cm*vF[,{i}]; w=w/w[2]; P=liftall(vK[{i}])')
            d=int(pari('poldegree(P,y)')) if 'y' in str(pari('variables(P)')) else 1
            coords=[]; trs=[]
            for p in PR:
                if d>1: coords.append([to24(tcoeffs(f'polcoef(liftall(w[{p}+1]),{j},y)',o),max(o,1)) for j in range(d)])
                else: coords.append([to24(tcoeffs(f'w[{p}+1]',o),max(o,1))])
                pari(f'a=w[{p}+1]; if(type(a)=="t_POLMOD" && variable(a.mod)==y, a=trace(a))')
                trs.append(to24(tcoeffs('a',o),max(o,1)))
            res['orbits'].append(dict(index=i,d=d,poly=str(pari('P')),coords=coords,traces=trs))
    pickle.dump(res,open(fn,'wb')); return res
if __name__=="__main__":
    t0=time.time()
    plan=pickle.load(open('plan.pkl','rb')) if os.path.exists('plan.pkl') else {}
    for bi,r in enumerate(reps):
        if bi not in plan:
            sp=[]
            for a in range(BD[8][r[0]][0]+1):
                for b in range(BD[9][r[1]][0]+1):
                    for c in range(BD[5][r[2]][0]+1):
                        M=2**a*3**b*5**c
                        if M<11: continue                      # no weight-2 newforms of level < 11 (Gamma_1 genus 0)
                        for conj in (False,True):
                            ks=charvec(M,r,conj)
                            if ks is not None and (M,tuple(ks)) not in sp: sp.append((M,tuple(ks)))
            plan[bi]=sp; pickle.dump(plan,open('plan.pkl','wb'))
    todo=sorted(set(x for sp in plan.values() for x in sp))
    left=[x for x in todo if not os.path.exists(f'spaces/{x[0]}_{"_".join(map(str,x[1]))}.pkl')]
    for (M,ks) in sorted(left):
        if time.time()-t0>230: break
        space(M,ks)
    left=[x for x in todo if not os.path.exists(f'spaces/{x[0]}_{"_".join(map(str,x[1]))}.pkl')]
    print("distinct candidate spaces:",len(todo),"; computed:",len(todo)-len(left),"; remaining:",len(left))
