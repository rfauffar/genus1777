# Stage 10: write the self-contained certificate (gzipped JSON) with all raw data needed by the standalone verifier.
import pickle, json, gzip, math, cypari2
from fractions import Fraction as Fr
from sympy import factorint
pari=cypari2.Pari(); pari.allocatemem(int(3*10**9)); pari('default(nbthreads,1)')
G=pickle.load(open('group.pkl','rb')); BL=pickle.load(open('blocks.pkl','rb')); R=pickle.load(open('traces.pkl','rb'))
plan=pickle.load(open('plan.pkl','rb')); Mt=pickle.load(open('match.pkl','rb')); EL=pickle.load(open('elliptic.pkl','rb'))
PR=R['PR']; ell=Mt['ell']; Z=Mt['Z']
fr=lambda x: str(Fr(x))
vec=lambda v: [fr(c) for c in v]
cert=dict(
  description="Genus-1777 curve X_Gamma, Gamma = {g in SL2(Z): g = lambda I mod 8 (lambda^2=1), g = +-(1 b;0 1) mod 3, g = +-I mod 5}",
  group=dict(G8=G['G8'],G9=G['G9'],G5=G['G5'],coset=G['coset']),
  cyclotomic=dict(L=24,basis="zeta_24^0..zeta_24^7 (reduced modulo Phi_24 = x^8 - x^4 + 1)"),
  primes=PR, ell=ell, zeta24_mod_ell=Z,
  blocks=[], spaces={}, orbits={})
used_keys=set()
for bi,b in enumerate(BL['blocks']):
    r=R['reps'][bi]; m=Mt['blocks'][bi]
    cert['blocks'].append(dict(index=bi,rep=list(r),members=[list(x) for x in b['members']],dim=str(b['dim']),
        candidate_spaces=[[M,list(ks)] for M,ks in plan[bi]],
        minus_N_p=[vec(tuple(-c for c in R['N'][p][bi])) for p in PR],
        multiplicities=[[k[0],list(k[1]),k[2],n] for k,n in m['mult'].items()],
        rank_witness_rows=m.get('witness',{}).get('rows',[]),det_mod_ell=m.get('witness',{}).get('det_mod_ell')))
    for k in m['mult']: used_keys.add(k)
for bi in plan:
    for M,ks in plan[bi]:
        key=f'{M}_{"_".join(map(str,ks))}'
        if key in cert['spaces']: continue
        s=pickle.load(open(f'spaces/{key}.pkl','rb'))
        # witness rows needed for this space: union over blocks using it
        rows=sorted(set(i for bj in plan if (M,ks) in plan[bj] for i in cert['blocks'][bj]['rank_witness_rows']))
        gens=[int(x) for x in pari(f'lift(znstar({M},1).gen)')]; cyc=[int(x) for x in pari(f'znstar({M},1).cyc')]
        cert['spaces'][key]=dict(M=M,ks=list(ks),gens=gens,cyc=cyc,char_order=s['o'],newspace_dim=s['dim'],
            orbits=[dict(index=o['index'],rel_degree=o['d'],field=o['poly'],
                         traces=[vec(t) for t in o['traces']],
                         coords_at_rows={str(i):[vec(c) for c in o['coords'][i]] for i in rows}) for o in s['orbits']])
# ellipticity data with coefficients up to the bounds used
def sturm(L): return math.ceil(L*math.prod(1+1/p for p in factorint(L))/6)
for k,v in EL.items():
    entry=dict(v); entry.pop('quer_pairs',None)
    entry['quer_pairs']=[list(x) for x in v.get('quer_pairs',[])]
    if v['type']!='rational':
        M=v['M']; ks=k.split('#')[0].split('_')[1:]; idx=int(k.split('#')[1])
        B=max([v.get('B',0)]+[t.get('checked_up_to',0) for t in v.get('twists',[])]+[2400])
        pari(f'G=znstar({M},1); chi={[int(x) for x in ks]}; mf=mfinit([{M},2,[G,chi]],0); S=mfsplit(mf); P=liftall(S[2][{idx}]); Cm=mfcoefs(mf,{B}); w=Cm*S[1][,{idx}]; w=w/w[2]; w=subst(liftall(w),t,-1)')
        entry['field_abs']=str(pari('subst(P,t,-1)'))
        entry['a_n']=[str(pari(f'w[{n}+1]')) for n in range(B+1)]
    entry['space']=k.split('#')[0]
    cert['orbits'][k]=entry
cert['gap_tables']={m:open(f'tables_{m}.txt').read() for m in (8,9,5)}
with gzip.open('cert1777_full.json.gz','wt') as f: json.dump(cert,f)
import os
print("certificate written:",os.path.getsize('cert1777_full.json.gz'),"bytes (gzipped);",len(cert['blocks']),"blocks,",len(cert['spaces']),"spaces,",len(cert['orbits']),"orbits")
