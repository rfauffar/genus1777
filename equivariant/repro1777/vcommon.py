# Shared pure-Python helpers for the verifiers: exact Z[zeta_24] arithmetic, GAP-table parsing, conjugacy classes.
import math, itertools
PHI24=[1,0,0,0,-1,0,0,0,1]                       # x^8 - x^4 + 1 (low degree first)
def red24(v):
    v=list(v)+[0]*max(0,48-len(v))
    for k in range(len(v)-1,7,-1):
        c=v[k]
        if c:
            for j in range(9): v[k-8+j]-=c*PHI24[j]
    return tuple(v[:8])
def from_coeffs(coeffs,L):
    v=[0]*24
    for k,c in enumerate(coeffs): v[(k*(24//L))%24]+=c
    return red24(v)
def mul24(a,b):
    v=[0]*16
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y: v[i+j]+=x*y
    return red24(v)
ONE=(1,0,0,0,0,0,0,0)
def mulm(x,y,n): a,b,c,d=x; e,f,g,h=y; return ((a*e+b*g)%n,(a*f+b*h)%n,(c*e+d*g)%n,(c*f+d*h)%n)
def invm(x,n): a,b,c,d=x; D=pow((a*d-b*c)%n,-1,n); return ((d*D)%n,(-b*D)%n,(-c*D)%n,(a*D)%n)
def load_tables(gap_text):
    out={}
    for m,txt in gap_text.items():
        m=int(m); reps=[]; sizes=[]; chars=[]
        for line in txt.splitlines():
            w=line.split()
            if not w: continue
            if w[0]=='REP': reps.append(tuple(int(t)%m for t in w[1:5])); sizes.append(int(w[5]))
            elif w[0]=='L': L=int(w[1])
            elif w[0]=='CHAR': chars.append([tuple(map(int,v.split(','))) for v in w[1:]])
        p=[q for q in (2,3,5,7) if m%q==0][0]
        G=[x for x in itertools.product(range(m),repeat=4) if (x[0]*x[3]-x[1]*x[2])%p]
        gens=[(1,1,0,1),(1,0,1,1)]+[(u,0,0,1) for u in range(1,m) if math.gcd(u,m)==1]
        cid={}; k=0
        for r in reps:                                   # classes = conjugacy orbits of the GAP representatives
            orb=[r]; cid[r]=k; i=0
            while i<len(orb):
                y=orb[i]; i+=1
                for g in gens:
                    z=mulm(mulm(g,y,m),invm(g,m),m)
                    if z not in cid: cid[z]=k; orb.append(z)
            assert len(orb)==sizes[k], "class size mismatch with GAP"
            k+=1
        assert len(cid)==len(G), "GAP classes do not cover the group"
        vals=[[from_coeffs(ch[c],L) for c in range(len(reps))] for ch in chars]
        out[m]=dict(cid=cid,sizes=sizes,vals=vals,order=len(G))
    return out
