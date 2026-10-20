from fractions import Fraction as F
from functools import lru_cache
from math import factorial,comb

@lru_cache(None)
def size(t):return 1+sum(size(x) for x in t)
@lru_cache(None)
def trees(n,ordered=True):
 if n==1:return ((),)
 out=[]
 def rec(rem,cs):
  if rem==0:out.append(tuple(cs));return
  for k in range(1,rem+1):
   for t in trees(k,ordered):
    if ordered or not cs or repr(t)>=repr(cs[-1]):rec(rem-k,cs+[t])
 rec(n-1,[])
 return tuple(out)
@lru_cache(None)
def forest_basis(n):
 out=[]
 def rec(rem,cs):
  if rem==0:out.append(tuple(cs));return
  for k in range(1,rem+1):
   for t in trees(k,False):
    if not cs or repr(t)>=repr(cs[-1]):rec(rem-k,cs+[t])
 rec(n,[])
 return out
@lru_cache(None)
def exact(t):
 # Owren alpha divided by factorial(|t|-1)
 k=0;a=1
 for child in t:
  k+=size(child);a*=comb(k-1,size(child)-1)*exact(child)*factorial(size(child)-1)
 return a/F(factorial(size(t)-1))

def cf_coeff(U,V,p):
 ts=[t for n in range(1,p+2) for t in trees(n,True)]
 g={t:F(t==()) for t in ts};K={t:F(0) for t in ts}
 for u,v in zip(U,V):
  K={t:u*K[t]+g[t] for t in ts}
  E={}
  for t in ts:
   z=v**len(t)/F(factorial(len(t)))
   for ch in t:z*=K[ch]
   E[t]=z
  g={t:sum(g[t[:r]]*E[t[r:]] for r in range(len(t)+1)) for t in ts}
 return g

def classic_coeff(U,V,p):
 n=len(U);q=[F(0)]*n;w=q.copy();A=[]
 for i in range(n):
  A.append(w.copy());q=[U[i]*x for x in q];q[i]+=1;w=[x+V[i]*y for x,y in zip(w,q)]
 b=w; ts=[t for k in range(1,p+1) for t in trees(k,False)];ph={}
 for t in ts:
  v=[F(1)]*n
  for child in t:
   av=[sum(A[i][j]*ph[child][j] for j in range(n)) for i in range(n)]
   v=[x*y for x,y in zip(v,av)]
  ph[t]=v
 return {t:sum(b[i]*ph[t][i] for i in range(n)) for t in ts}

if __name__=='__main__':
 print({p:len(trees(p+1,True)) for p in range(7)})
