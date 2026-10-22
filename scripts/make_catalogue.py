# Writes data/newform_orbits.csv: every Galois orbit of weight-2 newforms of level M | 4800
# with character 1 or chi_5, with its coefficient field, multiplicity sigma_0(4800/M) and
# the traces Tr(a_p) for p < 100.  Usage: python3 make_catalogue.py ../data/newform_orbits.csv
import cypari2, csv, json, sys
pari=cypari2.Pari(); pari.allocatemem(8*10**9)
out=sys.argv[1]
N=4800; divs=[d for d in range(1,N+1) if N%d==0]
sig0=lambda n: sum(1 for d in range(1,n+1) if n%d==0)
PR=[p for p in range(2,100) if pari.isprime(p)]
rows=[]
for M in divs:
    for chi in ([1,5] if M%5==0 else [1]):
        mf=pari(f'mfinit([{M},2,{chi}],0)')
        if int(pari.mfdim(mf))==0: continue
        for k,(F,P) in enumerate(zip(pari.mfeigenbasis(mf),pari.mffields(mf))):
            a=pari.mfcoefs(F,100); d=int(pari.poldegree(P))
            tr=[str(int(pari.trace(a[p])) if d>1 else int(a[p])) for p in PR]
            rows.append([M,'1' if chi==1 else 'chi5',k+1,d,str(pari.polredabs(P)) if d>1 else 'y',sig0(N//M)]+tr)
with open(out,'w',newline='') as fh:
    w=csv.writer(fh)
    w.writerow(['level','character','orbit_index','degree','field_polynomial','multiplicity_sigma0']+[f'trace_a{p}' for p in PR])
    w.writerows(rows)
print(len(rows), sum(int(r[3])*int(r[5]) for r in rows))
