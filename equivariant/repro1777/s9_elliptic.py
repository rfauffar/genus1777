# Stage 9: ellipticity certificate for every newform orbit used (exact arithmetic in the Hecke field via PARI polmods).
#  rational : trivial character and Hecke field Q.
#  CM by D  : a_n*(1-chi_D(n)) = 0 for all n <= B, gcd(n, M*D)=1, B = Sturm bound of level L'=lcm(M, D^2).
#  inner twists (F=Q): for every automorphism s of the Hecke field K_f there is chi_s with s(a_n) = chi_s(n) a_n.
#     * complex conjugation is automatically an inner twist with chi = conj(eps) (Petersson; no check needed);
#     * every other s is checked exactly for all n <= B, gcd(n, M*D)=1, B = Sturm bound of level lcm(M, D^2), chi_s = (D/.).
#     If the inner twists exhaust Aut(K_f/Q) (|Gamma| = [K_f:Q]), then F = Q (Ribet).
#  Brauer class (Quer): B_f ~ prod_k (t_k, D_k), t_k generators of the square classes of a_p^2/eps(p), D_k the discriminant whose
#     Kronecker pattern matches; Hilbert symbols recorded at 2, all primes dividing t_k D_k, and infinity.
import pickle, math, cypari2
from sympy import primerange, factorint
pari=cypari2.Pari(); pari.allocatemem(int(3*10**9)); pari('default(nbthreads,1)')
used=pickle.load(open('used_orbits.pkl','rb'))
def sturm(L): return math.ceil(L*math.prod(1+1/p for p in factorint(L))/6)
DISCS=[D for D in range(-480,481) if D not in (0,1) and int(pari(f'isfundamental({D})')) and set(factorint(abs(D)))<={2,3,5}]
def is0(expr): return int(pari(f'({expr})==0'))==1
rep={}
for (M,ks,idx),uses in used:
    key=f'{M}_{"_".join(map(str,ks))}#{idx}'
    pari(f'G=znstar({M},1); chi={list(ks)}; mf=mfinit([{M},2,[G,chi]],0); o=charorder(G,chi); S=mfsplit(mf); P=liftall(S[2][{idx}])')
    o=int(pari('o')); d=int(pari('poldegree(P,y)')) if 'y' in str(pari('variables(P)')) else 1
    if o<=1 and d==1: rep[key]=dict(type='rational',M=M); continue
    assert o<=2, "all characters here are real"
    Pabs=pari('subst(P,t,-1)') if o==2 else pari('P')        # Q(chi)=Q, K_f = Q[y]/Pabs
    def load(B):
        pari(f'Cm=mfcoefs(mf,{B}); w=Cm*S[1][,{idx}]; w=w/w[2]; w=subst(liftall(w),t,-1); Pa={Pabs}')
        return B
    Bmax=load(400)
    a=lambda n: f'Mod(w[{n}+1],Pa)' if d>1 else f'w[{n}+1]'
    # --- CM ---
    cm=None
    for D in [D for D in DISCS if D<0]:
        if all(is0(a(p)) for p in primerange(7,200) if M%p and int(pari(f'kronecker({D},{p})'))==-1):
            B=sturm(math.lcm(M,D*D))
            if B>Bmax: Bmax=load(B)
            if all(is0(f'{a(n)}*(1-kronecker({D},{n}))') for n in range(1,B+1) if math.gcd(n,M*abs(D))==1):
                cm=dict(D=D,B=B); break
    if cm: rep[key]=dict(type='CM',M=M,char_order=o,degree=d,**cm); continue
    # --- inner twists ---
    autos=[str(x) for x in pari(f'nfgaloisconj({Pabs})')] if d>1 else ['y']
    is_galois=len(autos)==d
    twists=[]
    real=int(pari(f'polsturm({Pabs})'))==d
    for s in autos:
        if s=='y': twists.append(dict(sigma=s,chi='trivial')); continue
        # complex conjugation (if K_f is CM-type / imaginary) is automatic when it is an automorphism: test via embedding
        found=None
        for D in [1]+DISCS:
            ok=True
            for p in primerange(7,120):
                if M%p==0 or (D!=1 and abs(D)%p==0): continue
                if not is0(f'Mod(subst(lift({a(p)}),y,{s})-kronecker({D},{p})*lift({a(p)}),Pa)'): ok=False; break
            if not ok: continue
            B=sturm(math.lcm(M,D*D)) if D!=1 else sturm(M)
            if B>Bmax: Bmax=load(B)
            ok=all(is0(f'Mod(subst(lift({a(n)}),y,{s})-kronecker({D},{n})*lift({a(n)}),Pa)') for n in range(1,B+1) if math.gcd(n,M*abs(D))==1)
            if ok: found=dict(sigma=s,chi=f'kronecker({D},.)',D=D,checked_up_to=B); break
        twists.append(found if found else dict(sigma=s,chi=None))
    if not (is_galois and all(t['chi'] is not None for t in twists)):
        rep[key]=dict(type='NOT_ELLIPTIC',M=M,char_order=o,degree=d,galois=is_galois,twists=twists); continue
    # --- Brauer class (Quer) ---
    sq={}
    for p in primerange(7,min(Bmax,2400)+1):
        if M%p==0: continue
        e=int(pari(f'chareval(G,chi,{p})*2'))%2 if o==2 else 0      # eps(p) = (-1)^e
        b2=pari(f'simplify(lift(Mod(lift({a(p)})^2,Pa)))') if d>1 else pari(f'{a(p)}^2')
        assert str(pari(f'type({b2})')) in ('t_INT','t_FRAC'), "a_p^2 not rational"
        b2=int(b2)*(-1)**e
        if b2: sq[p]=int(pari(f'core({b2})'))
    gens=[]; span={1:0}
    for v in sorted(set(sq.values()),key=abs):
        if v not in span:
            k=len(gens); gens.append(v); span.update({int(pari(f'core({s_*v})')):(m|(1<<k)) for s_,m in list(span.items())})
    tw=[]
    for k,t in enumerate(gens):
        Dk=next((D for D in DISCS if all(((span[sq[p]]>>k)&1)==(1 if int(pari(f'kronecker({D},{p})'))==-1 else 0) for p in sq)),None)
        tw.append((t,Dk))
    places=sorted({2}|{int(q) for t,D in tw if D is not None for q in factorint(2*abs(t)*abs(D))})
    symbols={str(q):[int(pari(f'hilbert({t},{D},{q})')) for t,D in tw] for q in places+[0]} if all(D is not None for t,D in tw) else None
    split=symbols is not None and all(math.prod(v)==1 for v in symbols.values())
    rep[key]=dict(type='inner-twist, F=Q',M=M,char_order=o,degree=d,field=str(Pabs),twists=twists,
                  quer_pairs=tw,hilbert_symbols=symbols,brauer_split=split)
pickle.dump(rep,open('elliptic.pkl','wb'))
from collections import Counter
print(Counter(v['type'] for v in rep.values()))
for k,v in rep.items():
    if v['type']!='rational': print(k,v['type'],v.get('D',''),[t.get('D') for t in v.get('twists',[])],v.get('quer_pairs',''),v.get('brauer_split',''))
