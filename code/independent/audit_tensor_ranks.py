"""Independent shuffle-rank audit via choosing interleaving positions."""


def require(condition, detail="verification obligation failed"):
    """A proof obligation that remains active under python -O."""
    if not condition:
        raise ValueError(detail)

from itertools import combinations,permutations
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json,math

def interleavings(a,b):
    counts=Counter();n=len(a)+len(b)
    for places in combinations(range(n),len(a)):
        ap=set(places);i=j=0;word=[]
        for k in range(n):
            if k in ap:word.append(a[i]);i+=1
            else:word.append(b[j]);j+=1
        counts[''.join(word)]+=1
    return counts

def rref(rows):
    M=[list(map(F,r)) for r in rows];pivs=[];k=0
    for j in range(len(M[0])):
        ii=next((i for i in range(k,len(M)) if M[i][j]),None)
        if ii is None:continue
        M[ii],M[k]=M[k],M[ii];z=M[k][j];M[k]=[x/z for x in M[k]]
        for i in range(len(M)):
            if i!=k:
                z=M[i][j];M[i]=[a-z*b for a,b in zip(M[i],M[k])]
        pivs.append(j);k+=1
        if k==len(M):break
    return M,pivs

def determinant(M):
    A=[r[:] for r in M];n=len(A);d=1;sgn=1
    for k in range(n-1):
        r=next(i for i in range(k,n) if A[i][k])
        if r!=k:A[k],A[r]=A[r],A[k];sgn=-sgn
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=A[i][j]*A[k][k]-A[i][k]*A[k][j]
                require((val%d==0), 'audit_tensor_ranks.py:39: check failed');A[i][j]=val//d
        d=A[k][k]
        for i in range(k+1,n):A[i][k]=0
    return sgn*A[-1][-1]

p=Path(__file__).resolve().parents[1]
data=json.loads((p/'tensor_ranks.json').read_text());out=[]
for obj in data:
    pattern=''.join(obj['multiplicity_pattern']);words=sorted(set(''.join(x) for x in permutations(pattern)))
    rows={f'P({a})':[int(w[0]==a) for w in words] for a in sorted(set(pattern))}
    for word in words:
        for i in range(1,len(word)):
            a,b=word[:i],word[i:];counts=interleavings(a,b)
            rows[f'S({a},{b})']=[counts[w] for w in words]
    selected=[rows[x] for x in obj['equation_labels']]
    require((words==obj['words']), 'audit_tensor_ranks.py:54: check failed');require((selected==obj['matrix']), 'audit_tensor_ranks.py:54: check failed')
    matrix,piv=rref(list(rows.values()));rank=len(piv)
    require((rank==obj['rank']), 'audit_tensor_ranks.py:56: check failed')
    rec={'pattern':pattern,'unknowns':len(words),'full_rank':rank,'all_equations':len(rows)}
    if rank==len(words):
        det=determinant(selected);require((str(det)==obj['determinant']), 'audit_tensor_ranks.py:59: check failed');rec['determinant']=det
    else:
        free=next(i for i in range(len(words)) if i not in piv);v=[F(0)]*len(words);v[free]=1
        for i,j in enumerate(piv):v[j]=-matrix[i][free]
        lcm=math.lcm(*(x.denominator for x in v));v=[int(x*lcm) for x in v]
        require((any(v) and all(sum(a*b for a,b in zip(row,v))==0 for row in rows.values())), 'audit_tensor_ranks.py:64: check failed')
        rec['null_vector']=v;rec['word_order']=words
    out.append(rec)
report={'status':'PASS','method':'Independent positional shuffle enumeration and Bareiss determinant','patterns':out}
print(json.dumps(report,indent=2))
