"""Independent audit of the Williamson degree-six coefficient certificate.

No candidate_optimisation Python module is imported. The JSON coefficient
certificate is data; method coefficients are evaluated independently by the
closed product formula and by Owren's direct exponential recursion.
"""


def require(condition, detail="verification obligation failed"):
    """A proof obligation that remains active under python -O."""
    if not condition:
        raise ValueError(detail)

from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
from itertools import combinations_with_replacement
import ast, json, random, hashlib

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'order6_symbolic_certificate.json'
@lru_cache(None)
def nodes(t):return 1+sum(nodes(c) for c in t)
@lru_cache(None)
def flow(t):
    z=F(1);total=0
    for c in t:
        total+=nodes(c);z*=flow(c)/total
    return z
@lru_cache(None)
def ordinary_factorial(t):return nodes(t)*prod(ordinary_factorial(c) for c in t)
def prod(values):
    z=F(1)
    for v in values:z*=v
    return z
@lru_cache(None)
def ordered(n):
    if n==1:return ((),)
    def forest(weight):
        if weight==0:yield ();return
        for first_weight in range(1,weight+1):
            for first in ordered(first_weight):
                for rest in forest(weight-first_weight):yield (first,)+rest
    return tuple(forest(n-1))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def mv(A,v):return [dot(row,v) for row in A]
def tableau(U,V):
    n=len(U)
    G=[[V[i]*prod(U[k] for k in range(j+1,i+1)) if j<=i else F(0) for j in range(n)] for i in range(n)]
    A=[[sum((G[k][j] for k in range(i)),F(0)) for j in range(n)] for i in range(n)]
    return A,G

def tensor(G,args):
    # Direct convolution of the chronological row exponential word series.
    states=[F(1)]+[F(0)]*len(args)
    for row in G:
        out=[]
        for k in range(len(args)+1):
            out.append(sum((states[j]*prod(dot(row,arg) for arg in args[j:k])/factorial(k-j) for j in range(k+1)),F(0)))
        states=out
    return states[-1]

def main():
    data=json.loads(CERT.read_text())
    records=[]
    for r in data['coefficients']:
        t=ast.literal_eval(r['tree'])
        terms=[(ast.literal_eval(z['forest']),F(z['coefficient'])) for z in r['polynomial']]
        records.append((r['order'],t,terms))
    want=set(t for n in range(1,8) for t in ordered(n))
    require((len(records)==len(want)==197), 'audit_order6.py:65: check failed')
    require((set(t for _,t,_ in records)==want), 'audit_order6.py:66: check failed')
    require((all(order==nodes(t)-1 for order,t,_ in records)), 'audit_order6.py:67: check failed')
    require((all(sum(nodes(x) for x in f)==order for order,t,p in records for f,_ in p)), 'audit_order6.py:68: check failed')
    require((all(sum(c*prod(1/ordinary_factorial(x) for x in f) for f,c in p)==flow(t) for _,t,p in records)), 'audit_order6.py:69: check failed')
    random.seed(314159)
    cases=[]
    for n in [1,2,3,5,8,12]:
        cases.append(([F(random.randint(-3,3),2) for _ in range(n)],[F(random.randint(-4,4),3) for _ in range(n)]))
    cases += [([F(0)]*5,[F(0)]*5),([F(0),F(-1),F(-1),F(-11),F(1,10)],[F(1,2),F(2,3),F(-1,2),F(-1,10),F(1,6)])]
    all_ts=sorted(want,key=lambda t:(nodes(t),repr(t)))
    prefix_checks=0;row_checks=0;tensor_checks=0
    required=set(x for _,_,p in records for f,_ in p for x in f)
    # Add every ordinary subtree needed to compute stage elementary weights.
    todo=list(required)
    while todo:
        for c in todo.pop():
            if c not in required:required.add(c);todo.append(c)
    for U,V in cases:
        n=len(U);A,G=tableau(U,V)
        phi={}
        for t in sorted(required,key=lambda x:(nodes(x),repr(x))):
            ph=[F(1)]*n
            for c in t:
                ac=mv(A,phi[c]);ph=[x*y for x,y in zip(ph,ac)]
            phi[t]=ph
        b=[F(0)]*n
        state={t:F(t==()) for t in all_ts};delta={t:F(0) for t in all_ts}
        for stage in range(n+1):
            moments={t:dot(b,phi[t]) for t in required}
            for _,t,p in records:
                predicted=sum((c*prod(moments[x] for x in f) for f,c in p),F(0))
                require((predicted==state[t]), (n,stage,t,predicted,state[t]))
                prefix_checks+=1
            if stage==n:break
            delta={t:U[stage]*delta[t]+state[t] for t in all_ts}
            exp_coeff={t:V[stage]**len(t)*prod(delta[c] for c in t)/factorial(len(t)) for t in all_ts}
            state={t:sum((state[t[:j]]*exp_coeff[t[j:]] for j in range(len(t)+1)),F(0)) for t in all_ts}
            b=[x+y for x,y in zip(b,G[stage])]
        for _ in range(3):
            u=[F(random.randint(-4,4),5) for _ in range(n)]
            v=[F(random.randint(-4,4),7) for _ in range(n)]
            av=mv(A,v);bv=dot(b,v)
            R=[[ (bv-av[i] if i==j else F(0))+v[i]*A[i][j] for j in range(n)] for i in range(n)]
            for i,row in enumerate(G):
                lhs=[sum((row[k]*R[k][j] for k in range(n)),F(0)) for j in range(n)]
                rhs=[(bv-av[i])*row[j]+dot(row,v)*A[i][j] for j in range(n)]
                require((lhs==rhs), 'audit_order6.py:112: check failed');row_checks+=1
            ru=u
            for m in range(1,8):
                require((tensor(G,[u]+[v]*(m-1))==dot(b,ru)/factorial(m)), 'audit_order6.py:115: check failed')
                tensor_checks+=1;ru=mv(R,ru)
    report={'status':'PASS','certificate_sha256':hashlib.sha256(CERT.read_bytes()).hexdigest(),
            'ordered_coefficients':len(records),'coverage_by_order':[len(ordered(n)) for n in range(1,8)],
            'homogeneity_checks':sum(len(p) for _,_,p in records),'independent_exact_flow_checks':len(records),
            'coefficient_lists':len(cases),'prefix_coefficient_checks':prefix_checks,'general_vector_row_checks':row_checks,
            'general_vector_tensor_checks':tensor_checks,'all_arithmetic':'fractions.Fraction; no imported candidate code'}
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
