# Stage 6: block representatives, central characters, and EXACT Frobenius traces N_p(rho) in Q(zeta_24) for 7<=p<=2400, p !| 360,
# via the Duke-Toth / Eichler-Selberg class-number formula:
#   N_p(rho) = sum_{t^2<4p} sum_{b^2 | t^2-4p, D=(t^2-4p)/b^2 = 0,1 mod 4} h(D)/w(D) chi_rho(A_{t,b})
#            + 1/(2N) sum_{k mod N} sum_{e=+-1} chi_rho( e*(p, k; 0, 1) ),
#   A_{t,b} = ((t+b*dl)/2, b; b*(D-dl)/4, (t-b*dl)/2), dl = D mod 2 (trace t, determinant p).
# The trace of Frobenius on the block's Hecke module is -N_p(rho).
import pickle, math, os, time, sys
from fractions import Fraction
import cypari2, cyclo
from sympy import primerange
pari=cypari2.Pari()
C=pickle.load(open('classes.pkl','rb')); T=pickle.load(open('tables.pkl','rb')); B=pickle.load(open('blocks.pkl','rb')); BD=pickle.load(open('bounds.pkl','rb'))
N=360; PR=[p for p in primerange(7,2401) if N%p]
chi={m:[[cyclo.from_coeffs(T[m]['ex'][i][k],T[m]['L']) for k in range(len(C[m][1]))] for i in range(len(T[m]['ex']))] for m in (8,9,5)}
def lvl(r): return 2**BD[8][r[0]][0]*3**BD[9][r[1]][0]*5**BD[5][r[2]][0]
reps=[min(b['members'],key=lambda r:(lvl(r),r)) for b in B['blocks']]
def cls3(x): return (C[8][0][tuple(v%8 for v in x)],C[9][0][tuple(v%9 for v in x)],C[5][0][tuple(v%5 for v in x)])
def hw(D): return int(pari.qfbclassno(D)),(6 if D==-3 else (4 if D==-4 else 2))
fn='traces.pkl'; R=pickle.load(open(fn,'rb')) if os.path.exists(fn) else {'reps':reps,'PR':PR,'N':{}}
t0=time.time(); cusp_cache={}
for p in PR:
    if p in R['N']: continue
    if time.time()-t0>240: break
    W={}
    Tb=2*math.isqrt(p)+1
    for t in range(-Tb,Tb+1):
        disc=t*t-4*p
        if disc>=0: continue
        b=1
        while b*b<=-disc:
            if disc%(b*b)==0:
                D=disc//(b*b)
                if D%4 in (0,1):
                    h,w=hw(D); dl=D%2
                    k=cls3(((t+b*dl)//2,b,b*(D-dl)//4,(t-b*dl)//2)); W[k]=W.get(k,0)+Fraction(h,w)
            b+=1
    key=p%N
    if key not in cusp_cache:
        Wc={}
        for kk in range(N):
            for e in (1,-1):
                k=cls3((e*p,e*kk,0,e)); Wc[k]=Wc.get(k,0)+Fraction(1,2*N)
        cusp_cache[key]=Wc
    for k,v in cusp_cache[key].items(): W[k]=W.get(k,0)+v
    vals=[]
    for r in reps:
        s=[Fraction(0)]*cyclo.DEG
        for (a,b2,c),wt in W.items():
            v=cyclo.mul(cyclo.mul(chi[8][r[0]][a],chi[9][r[1]][b2]),chi[5][r[2]][c])
            for j in range(cyclo.DEG): s[j]+=wt*v[j]
        vals.append(tuple(s))
    R['N'][p]=vals
pickle.dump(R,open(fn,'wb'))
print("blocks:",len(reps),"; primes done:",len(R['N']),"of",len(PR))
