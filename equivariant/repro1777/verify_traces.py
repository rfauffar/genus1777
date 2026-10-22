#!/usr/bin/env python3
# V6: recompute the exact block traces -N_p(rho) from the GAP character tables and class numbers computed here
#     (reduced primitive binary quadratic forms), and compare with the certificate.  Pure Python.
import gzip, json, sys, math, time
from fractions import Fraction as Fr
from vcommon import *
cert=json.load(gzip.open(sys.argv[1] if len(sys.argv)>1 else 'cert1777_full.json.gz','rt'))
TB=load_tables(cert['gap_tables']); N=360
_h={}
def classno(D):
    # number of SL2(Z)-classes of primitive positive definite forms of discriminant D<0
    if D in _h: return _h[D]
    h=0; a=1
    while 3*a*a<=-D:
        for b in range(-a+1,a+1):
            if (b*b-D)%(4*a): continue
            c=(b*b-D)//(4*a)
            if c<a: continue
            if b<0 and c==a: continue
            if math.gcd(math.gcd(a,abs(b)),c)!=1: continue
            h+=1
        a+=1
    _h[D]=h; return h
def cls3(x): return tuple(TB[m]['cid'][tuple(v%m for v in x)] for m in (8,9,5))
def chi(r,k): return mul24(mul24(TB[8]['vals'][r[0]][k[0]],TB[9]['vals'][r[1]][k[1]]),TB[5]['vals'][r[2]][k[2]])
reps=[tuple(b['rep']) for b in cert['blocks']]; PR=cert['primes']; bad=[]; cusp={}
start=int(sys.argv[2]) if len(sys.argv)>2 else 0; t0=time.time(); done=0
for ip in range(start,len(PR)):
    p=PR[ip]
    if time.time()-t0>250: print(f"(time limit: stopped before prime index {ip}; rerun with start index {ip})"); break
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
                    w=6 if D==-3 else (4 if D==-4 else 2); dl=D%2
                    k=cls3(((t+b*dl)//2,b,b*(D-dl)//4,(t-b*dl)//2)); W[k]=W.get(k,0)+Fr(classno(D),w)
            b+=1
    if p%N not in cusp:
        Wc={}
        for kk in range(N):
            for e in (1,-1):
                k=cls3((e*p,e*kk,0,e)); Wc[k]=Wc.get(k,0)+Fr(1,2*N)
        cusp[p%N]=Wc
    for k,v in cusp[p%N].items(): W[k]=W.get(k,0)+v
    for bi,r in enumerate(reps):
        s=[Fr(0)]*8
        for k,wt in W.items():
            v=chi(r,k)
            for j in range(8): s[j]+=wt*v[j]
        if [-x for x in s]!=[Fr(x) for x in cert['blocks'][bi]['minus_N_p'][ip]]: bad.append((bi,p))
    done+=1
print(f"[{'PASS' if not bad else 'FAIL'}] V6 traces recomputed from the GAP tables for primes index {start}..{start+done-1} ({done} primes x {len(reps)} blocks) {bad[:3]}")
