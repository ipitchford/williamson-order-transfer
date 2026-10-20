"""Exact symbolic reduction of ordered coefficients using proven tensor identities.

All expressions are polynomials in ordinary rooted-tree moments M_tau.
No fitting or floating point enters the reduction.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict,Counter
from itertools import permutations,product
from math import factorial
import json
from tree_probe import trees,size,exact


def canon(cs):return tuple(sorted(cs,key=repr))
def clean(p):return {m:c for m,c in p.items() if c}
def add(p,q,scale=F(1)):
 r=dict(p)
 for m,c in q.items():r[m]=r.get(m,F(0))+scale*c
 return clean(r)
def scale(p,c):return clean({m:c*x for m,x in p.items()})
def mul(p,q):
 r=defaultdict(F)
 for m,c in p.items():
  for n,d in q.items():r[canon(m+n)]+=c*d
 return clean(r)
def one():return {():F(1)}

def r_op(expr,v):
 # expr keys: (moment forest, ordinary stage tree).
 out=defaultdict(F)
 for (f,x),c in expr.items():
  out[(canon(f+(v,)),x)]+=c
  out[(f,canon(x+(v,))) ]-=c # -(A Phi_v) * Phi_x
  out[(f,canon(v+(x,))) ]+=c # Phi_v * (A Phi_x)
 return clean(out)
def b_contract(expr):
 out=defaultdict(F)
 for (f,x),c in expr.items():out[canon(f+(x,))]+=c
 return clean(out)

@lru_cache(None)
def perms(t):return tuple(sorted(set(permutations(t)),key=repr))
@lru_cache(None)
def shuffles(a,b):
 if not a:return {b:1}
 if not b:return {a:1}
 out=defaultdict(int)
 for t,k in shuffles(a[1:],b).items():out[(a[0],)+t]+=k
 for t,k in shuffles(a,b[1:]).items():out[(b[0],)+t]+=k
 return dict(out)

FAMILY_LOG=[]
FAMILY_DATA={}
@lru_cache(None)
def tensor_family(ms):
 k=len(ms);ps=perms(ms);ix={p:i for i,p in enumerate(ps)};eqs=[]
 if k==1:
  answer={ms:{(ms[0],):F(1)}};FAMILY_DATA[ms]=answer;return answer
 # Polarised form of t_k(u,v,...,v)=b R_v^(k-1)u/k!.
 for first in sorted(set(ms),key=repr):
  tail=list(ms);tail.remove(first);tps=perms(tuple(tail));row=[F(0)]*len(ps);rhs={}
  for tp in tps:
   row[ix[(first,)+tp]]+=1
   expr={((),first):F(1)}
   for v in tp:expr=r_op(expr,v)
   rhs=add(rhs,b_contract(expr),F(1,factorial(k)))
  eqs.append((row,rhs))
 # Chen's shuffle product; all proper cuts give small exact linear systems.
 for p in ps:
  for j in range(1,k):
   a=p[:j];b=p[j:];row=[F(0)]*len(ps)
   for q,c in shuffles(a,b).items():row[ix[q]]+=c
   rhs=mul(tensor(a),tensor(b));eqs.append((row,rhs))
 piv={}
 for row,rhs in eqs:
  for i in sorted(piv):
   if row[i]:
    c=row[i];pr,pb=piv[i];row=[x-c*y for x,y in zip(row,pr)];rhs=add(rhs,pb,-c)
  i=next((i for i,x in enumerate(row) if x),None)
  if i is None:
   # A residual here is an identity among classical moments specific to low storage.
   # Store it; do not silently assume generic independence of those moments.
   if rhs:RELATIONS.append(rhs)
   continue
  c=row[i];row=[x/c for x in row];rhs=scale(rhs,1/c);piv[i]=(row,rhs)
 if len(piv)<len(ps):raise RuntimeError(f'Unresolved tensor family weight {sum(size(t) for t in ms)}: {ms}, rank {len(piv)}/{len(ps)}')
 for i in sorted(piv,reverse=True):
  row,rhs=piv[i]
  for j in range(i+1,len(ps)):
   if row[j]:rhs=add(rhs,piv[j][1],-row[j]);row[j]=F(0)
  piv[i]=(row,rhs)
 FAMILY_LOG.append({'weight':sum(size(t) for t in ms),'length':k,'permutations':len(ps),'rank':len(piv)})
 answer={p:piv[i][1] for i,p in enumerate(ps)}
 FAMILY_DATA[ms]=answer
 return answer
RELATIONS=[]
def tensor(args):return tensor_family(canon(args))[args]

@lru_cache(None)
def gpoly(t):
 if not t:return one()
 choices=[]
 for child in t:
  choices.append([(canon(f),c) for f,c in gpoly(child).items()])
 out={}
 for combo in product(*choices):
  args=tuple(x[0] for x in combo);c=F(1)
  for _,d in combo:c*=d
  out=add(out,tensor(args),c)
 return out

@lru_cache(None)
def tree_factorial(t):
 z=size(t)
 for ch in t:z*=tree_factorial(ch)
 return z

def exact_moment_eval(poly):
 out=F(0)
 for f,c in poly.items():
  for t in f:c/=tree_factorial(t)
  out+=c
 return out

def serial(poly):return [{'forest':repr(f),'coefficient':str(c)} for f,c in sorted(poly.items(),key=lambda x:repr(x[0]))]

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,default=6);parser.add_argument('--output',default=None);args=parser.parse_args()
 records=[];by_order={};fail=None
 for p in range(args.order+1):
  n=0;bad=[]
  try:
   for t in trees(p+1,True):
    poly=gpoly(t);value=exact_moment_eval(poly)
    if value!=exact(t):bad.append({'tree':repr(t),'got':str(value),'wanted':str(exact(t))})
    records.append({'order':p,'tree':repr(t),'polynomial':serial(poly),'exact_moment_value':str(value),'exact_flow_value':str(exact(t))});n+=1
  except RuntimeError as e:fail=str(e)
  by_order[p]={'ordered_trees':n,'errors':bad};print('order',p,'trees',n,'errors',len(bad),flush=True)
  if fail:print(fail,flush=True);break
 report={'status':'PASS' if not fail and not any(x['errors'] for x in by_order.values()) else 'STOP','order_requested':args.order,'by_order':by_order,
 'method':'Exact polarised repeated-tail tensor formula plus shuffle identities; no numerical fitting',
 'tensor_families':len(FAMILY_LOG),'maximum_family_size':max(x['permutations'] for x in FAMILY_LOG),'relations_encountered':len(RELATIONS),
 'relations_at_exact_moments':[str(exact_moment_eval(p)) for p in RELATIONS], 'unresolved':fail,'coefficients':records,
 'tensor_formulas':[{'multiset':repr(ms),'weight':sum(size(t) for t in ms),
 'formulas':[{'word':repr(w),'polynomial':serial(p)} for w,p in sorted(data.items(),key=lambda x:repr(x[0]))]}
 for ms,data in sorted(FAMILY_DATA.items(),key=lambda x:(sum(size(t) for t in x[0]),repr(x[0])))]}
 if args.output:
  with open(args.output,'w') as f:json.dump(report,f,indent=2)
 print(json.dumps({k:v for k,v in report.items() if k not in ['coefficients','relations_at_exact_moments','tensor_formulas']},indent=2))
