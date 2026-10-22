# Stage 4: group the contributing representations into determinant-twist classes (blocks); block dimensions.
import pickle, math
from fractions import Fraction
import cyclo
C=pickle.load(open('classes.pkl','rb')); T=pickle.load(open('tables.pkl','rb')); G=pickle.load(open('group.pkl','rb'))
contrib=G['contrib']; phiN=G['phiN']
def is_det(m,i):
    cid=C[m][0]
    if cyclo.from_coeffs(T[m]['ex'][i][cid[(1,0,0,1)]],T[m]['L'])!=cyclo.ONE: return False
    return all(cyclo.from_coeffs(T[m]['ex'][i][cid[x]],T[m]['L'])==cyclo.ONE for x in cid if (x[0]*x[3]-x[1]*x[2])%m==1)
def twist_table(m):
    # index of rho (x) delta for every irreducible rho and determinant character delta (exact comparison of characters)
    n=len(T[m]['ex']); cid=C[m][0]; L=T[m]['L']
    vals=[[cyclo.from_coeffs(T[m]['ex'][i][k],L) for k in range(len(C[m][1]))] for i in range(n)]
    dets=[i for i in range(n) if is_det(m,i)]
    key={tuple(v):i for i,v in enumerate(vals)}
    tw={}
    for i in range(n):
        for d in dets: tw[(i,d)]=key[tuple(cyclo.mul(a,b) for a,b in zip(vals[i],vals[d]))]
    return dets,tw
TW={m:twist_table(m) for m in (8,9,5)}
def orbit(r):
    return sorted(set((TW[8][1][(r[0],a)],TW[9][1][(r[1],b)],TW[5][1][(r[2],c)]) for a in TW[8][0] for b in TW[9][0] for c in TW[5][0]))
blocks={}
for r in contrib:
    o=tuple(orbit(r)); blocks.setdefault(o[0],set()).add(r)
bl=[]
for c,mem in sorted(blocks.items()):
    d=Fraction(sum(contrib[r] for r in mem),phiN); bl.append(dict(canon=c,members=sorted(mem),dim=d))
print("blocks (determinant-twist classes):",len(bl),"; sum of block dimensions:",sum(b['dim'] for b in bl))
print("block dimensions:",sorted(int(b['dim']) for b in bl))
pickle.dump(dict(blocks=bl,TW=TW),open('blocks.pkl','wb'))
