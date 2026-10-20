from itertools import permutations
from fractions import Fraction as F
from functools import lru_cache
import json

@lru_cache(None)
def shuffles(a,b):
 if not a:return {b:1}
 if not b:return {a:1}
 out={}
 for t,c in shuffles(a[1:],b).items():out[(a[0],)+t]=out.get((a[0],)+t,0)+c
 for t,c in shuffles(a,b[1:]).items():out[(b[0],)+t]=out.get((b[0],)+t,0)+c
 return out

def certificate(ms):
 ps=sorted(set(permutations(ms)));ix={p:i for i,p in enumerate(ps)};rows=[];labels=[]
 for first in sorted(set(ms)):
  rows.append([int(p[0]==first) for p in ps]);labels.append('P('+first+')')
 for p in ps:
  for cut in range(1,len(ms)):
   left,right=p[:cut],p[cut:];sh=shuffles(left,right)
   rows.append([sh.get(q,0) for q in ps]);labels.append('S('+''.join(left)+','+''.join(right)+')')
 piv={};selected=[];chosen=[]
 for row,label in zip(rows,labels):
  v=list(map(F,row))
  for i in sorted(piv):
   z=v[i]
   if z:v=[a-z*b for a,b in zip(v,piv[i])]
  i=next((i for i,z in enumerate(v) if z),None)
  if i is not None:
   z=v[i];v=[a/z for a in v];piv[i]=v;selected.append(row);chosen.append(label)
 if len(piv)==len(ps):
  mat=[list(map(F,r)) for r in selected];det=F(1)
  for c in range(len(ps)):
   j=next(j for j in range(c,len(ps)) if mat[j][c])
   if j!=c:mat[c],mat[j]=mat[j],mat[c];det=-det
   z=mat[c][c];det*=z;mat[c]=[x/z for x in mat[c]]
   for j in range(c+1,len(ps)):
    z=mat[j][c];mat[j]=[a-z*b for a,b in zip(mat[j],mat[c])]
 else:det=None
 return {'multiplicity_pattern':ms,'words':[''.join(p) for p in ps],'rank':len(piv),'equation_labels':chosen,'matrix':selected,'determinant':str(det) if det is not None else None}
if __name__=='__main__':
 data=[certificate(tuple(ms)) for ms in ['ab','aab','abc','aaab','aabb','aaaab','aaabb']]
 print(json.dumps(data,indent=2))
