# Stage 2: conjugacy-class maps of GL2(Z/m) (m=8,9,5) matched to the GAP classes of Stage 1; exact character tables.
# Writes classes.pkl = {m: (cid: matrix->class index, sizes, |G|)}, tables.pkl = {m: {'L','ex' (exact coeff vectors),'num'}}.
import itertools, math, cmath, pickle, numpy as np
def mul(x,y,n): a,b,c,d=x; e,f,g,h=y; return ((a*e+b*g)%n,(a*f+b*h)%n,(c*e+d*g)%n,(c*f+d*h)%n)
def inv(x,n):
    a,b,c,d=x; D=pow((a*d-b*c)%n,-1,n); return ((d*D)%n,(-b*D)%n,(-c*D)%n,(a*D)%n)
def units(n): return [u for u in range(1,n) if math.gcd(u,n)==1]
C={}; T={}
for m in (8,9,5):
    p=[q for q in (2,3,5) if m%q==0][0]
    G=[x for x in itertools.product(range(m),repeat=4) if (x[0]*x[3]-x[1]*x[2])%p]
    gens=[(1,1,0,1),(1,0,1,1)]+[(u,0,0,1) for u in units(m)]          # generate GL2(Z/m)
    cid={}; sizes=[]
    for x in G:                                                          # conjugacy classes by orbit enumeration
        if x in cid: continue
        k=len(sizes); orb=[x]; cid[x]=k; i=0
        while i<len(orb):
            y=orb[i]; i+=1
            for g in gens:
                z=mul(mul(g,y,m),inv(g,m),m)
                if z not in cid: cid[z]=k; orb.append(z)
        sizes.append(len(orb))
    reps=[]; gsz=[]; chars=[]
    for line in open(f'tables_{m}.txt'):
        w=line.split()
        if w[0]=='REP': reps.append(tuple(int(t)%m for t in w[1:5])); gsz.append(int(w[5]))
        elif w[0]=='L': L=int(w[1])
        elif w[0]=='CHAR': chars.append([tuple(map(int,v.split(','))) for v in w[1:]])
    perm=[cid[r] for r in reps]
    assert sorted(perm)==list(range(len(sizes))), "class matching failed"
    assert all(sizes[perm[j]]==gsz[j] for j in range(len(reps))), "class sizes disagree with GAP"
    # re-index our classes to GAP's order
    newid={perm[j]:j for j in range(len(reps))}
    cid={x:newid[k] for x,k in cid.items()}; sizes=gsz
    z=[cmath.exp(2j*math.pi*k/L) for k in range(L)]
    num=np.array([[sum(c*z[k] for k,c in enumerate(v)) for v in ch] for ch in chars])
    orth=(num*np.array(sizes))@num.conj().T/len(G)
    assert np.allclose(orth,np.eye(len(chars))), "orthogonality fails"
    C[m]=(cid,sizes,len(G)); T[m]={'L':L,'ex':chars,'num':num}
    print(f"GL2(Z/{m}): order {len(G)}, {len(sizes)} classes (matched to GAP, sizes agree), orthogonality verified")
pickle.dump(C,open('classes.pkl','wb')); pickle.dump(T,open('tables.pkl','wb'))
