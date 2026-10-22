# Stage 5: local conductor bounds (corrected lemma): for n>=1, 0<=k<=e, 0<=n-k<=e, if rho has a vector fixed by
#   H_{k,n} = {(a, p^k b; p^(n-k) c, d) mod p^e : det unit, d = 1 mod p^min(n,e)}  then a(pi_p) <= n.
# n=0 (spherical): H = all of GL2(Z/p^e); only the trivial representation has invariants, giving a(pi)=0.
# Fallback without invariants: 2e (twist-minimal bound for K(p^e)-types).  Also records the central character.
import itertools, pickle
import cyclo
C=pickle.load(open('classes.pkl','rb')); T=pickle.load(open('tables.pkl','rb'))
def H(m,p,e,k,n):
    S=set()
    for a,b,c,d in itertools.product(range(m),repeat=4):
        if n>0 and (d-1)%(p**min(n,e)): continue
        x=(a,(p**k*b)%m,(p**(n-k)*c)%m if n-k<=e else 0,d)
        if (x[0]*x[3]-x[1]*x[2])%p: S.add(x)
    return S
bounds={}
for m,p,e in ((8,2,3),(9,3,2),(5,5,1)):
    cid=C[m][0]; L=T[m]['L']; nirr=len(T[m]['ex']); best=[None]*nirr; how=[None]*nirr
    for n in range(0,2*e+1):
        for k in range(0,min(n,e)+1):
            if n-k>e: continue
            Hs=H(m,p,e,k,n); cnt={}
            for x in Hs: cnt[cid[x]]=cnt.get(cid[x],0)+1
            for i in range(nirr):
                if best[i] is not None: continue
                s=cyclo.ZERO
                for cl,c in cnt.items():
                    v=cyclo.from_coeffs(T[m]['ex'][i][cl],L); s=cyclo.add(s,tuple(c*t for t in v))
                assert cyclo.is_rational(s) and s[0]%len(Hs)==0
                if s[0]>0: best[i]=n; how[i]=(k,n)
    bounds[m]=[(b if b is not None else 2*e, h) for b,h in zip(best,how)]
    from collections import Counter
    print(m,"conductor-bound histogram:",sorted(Counter(b for b,_ in bounds[m]).items()))
pickle.dump(bounds,open('bounds.pkl','wb'))
