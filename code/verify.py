"""Exact rational checks for the order-six transfer manuscript.

This complements the analytic induction; finite examples do not replace it.
Requires Python 3.10+ and NumPy. No network calls.
"""


def require(condition, detail="verification obligation failed"):
    """A proof obligation that remains active under python -O."""
    if not condition:
        raise ValueError(detail)

from fractions import Fraction as F
from math import factorial
import ast, json, random
from pathlib import Path
import numpy as np
from tree_probe import cf_coeff, classic_coeff, exact, size


def tableau(U,V):
    n=len(U);q=np.zeros(n,dtype=object);w=q.copy();a=[];G=[]
    for i in range(n):
        a.append(w.copy());q=U[i]*q;q[i]+=1;G.append(V[i]*q.copy());w=w+V[i]*q
    return np.array(a,dtype=object),w,np.array(G,dtype=object)


def tensor(G,args):
    # Dynamic product of exponentials, preserving chronological order.
    k=len(args);T=[F(1)]+[F(0)]*k
    for row in G:
        new=T.copy()
        for j in range(1,k+1):
            for r in range(1,j+1):
                z=T[j-r]/F(factorial(r))
                for u in args[j-r:j]:z*=row@u
                new[j]+=z
        T=new
    return T[-1]


def identity_checks():
    random.seed(20261019);lists=[]
    for n in [1,2,3,5,8,12]:
        for _ in range(3):
            U=[F(random.randint(-3,3),2) for _ in range(n)]
            V=[F(random.randint(-3,3),3) for _ in range(n)]
            lists.append((U,V))
    lists+=[([F(0)]*4,[F(0)]*4),([F(1)]*3,[F(1)]*3)]
    checks=0;general_checks=0
    for U,V in lists:
        A,b,G=tableau(U,V);n=len(U);one=np.ones(n,dtype=object);c=A@one;C=b@one;D=np.diag(c)
        R=C*np.eye(n,dtype=object)-D+A
        for i in range(n):
            require((list(G[i]@(D-A))==list(c[i]*G[i]-(G[i]@one)*A[i])), 'verify.py:48: check failed');checks+=1
        for _ in range(3):
            u=np.array([F(random.randint(-3,3),5) for _ in range(n)],dtype=object);v=u.copy()
            for m in range(1,9):
                require((tensor(G,[u]+[one]*(m-1))==(b@v)/F(factorial(m))), 'verify.py:52: check failed');checks+=1
                v=R@v
        for _ in range(3):
            u=np.array([F(random.randint(-3,3),5) for _ in range(n)],dtype=object)
            v=np.array([F(random.randint(-3,3),7) for _ in range(n)],dtype=object)
            Rv=(b@v)*np.eye(n,dtype=object)-np.diag(A@v)+np.diag(v)@A
            for i in range(n):
                require((G[i]@((A@v)*u-v*(A@u))==(A[i]@v)*(G[i]@u)-(G[i]@v)*(A[i]@u)), 'verify.py:59: check failed')
                general_checks+=1
            power=u.copy()
            for m in range(1,9):
                require((tensor(G,[u]+[v]*(m-1))==(b@power)/F(factorial(m))), 'verify.py:63: check failed')
                power=Rv@power;general_checks+=1
    return {'coefficient_lists':len(lists),'unit_tail_rational_checks':checks,
            'general_vector_rational_checks':general_checks,'maximum_tensor_degree':8}


def certificate_example_checks():
    records=json.loads((Path(__file__).resolve().parent/'order6_symbolic_certificate.json').read_text())['coefficients']
    rng=random.Random(641926);checks=0
    for n in [2,3,5,8]:
        U=[F(rng.randint(-3,3),2) for _ in range(n)]
        V=[F(rng.randint(-3,3),3) for _ in range(n)]
        moments=classic_coeff(U,V,6);ordered=cf_coeff(U,V,6)
        for record in records:
            value=F(0)
            for term in record['polynomial']:
                coefficient=F(term['coefficient'])
                for tree in ast.literal_eval(term['forest']):coefficient*=moments[tree]
                value+=coefficient
            require((value==ordered[ast.literal_eval(record['tree'])]), 'verify.py:82: check failed');checks+=1
    return {'arbitrary_rational_lists':4,'direct_ordered_series_comparisons':checks}


def p4_checks():
    U=list(map(F,[0,-1,-1,-11]))+[F(1,10)]
    V=[F(1,2),F(2,3),F(-1,2),F(-1,10),F(1,6)]
    A,b,G=tableau(U,V);one=np.ones(5,dtype=object);c=A@one;ac=A@c
    got=[b@one,b@c,b@(c*c),b@ac,b@(c**3),b@(c*ac),b@A@(c*c),b@A@ac]
    want=list(map(F,['1','1/2','1/3','1/6','1/4','1/8','1/12','1/24']))
    require((got==want), 'verify.py:92: check failed')
    ordered=[tensor(G,[c,one]),tensor(G,[c,one,one]),tensor(G,[c*c,one]),tensor(G,[ac,one]),sum(b[r]*tensor(G[:r],[c,one]) for r in range(5))]
    require((ordered==list(map(F,['1/6','1/24','1/12','1/24','1/24']))), 'verify.py:94: check failed')
    series=cf_coeff(U,V,5)
    require((all(z==exact(t) for t,z in series.items() if size(t)<=5)), 'verify.py:96: check failed')
    defects=sum(z!=exact(t) for t,z in series.items() if size(t)==6)
    return U,V,{'classical_conditions':[str(v) for v in got],'additional_conditions':[str(v) for v in ordered],
                'all_ordered_coefficients_through_degree_four':sum(size(t)<=5 for t in series),'fifth_degree_defects':defects,
                'A':[[str(x) for x in row] for row in A],'b':[str(x) for x in b],'c':[str(x) for x in c]}

# Direct polynomial matrix arithmetic, independent of ordered rooted trees.
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def zmat():return [[F(0),F(0)],[F(0),F(0)]]
def madd(A,B):return [[A[i][j]+B[i][j] for j in range(2)] for i in range(2)]
def mscale(A,q):return [[q*x for x in row] for row in A]
def pzero(p):return [zmat() for _ in range(p+1)]
def padd(A,B):return [madd(a,b) for a,b in zip(A,B)]
def pscale(A,q):return [mscale(a,q) for a in A]
def pmul(A,B):
 p=len(A)-1;out=pzero(p)
 for n in range(p+1):
  for i in range(n+1):out[n]=madd(out[n],mm(A[i],B[n-i]))
 return out
I=[[F(1),F(0)],[F(0),F(1)]]
X=[[F(0),F(1)],[F(-1),F(0)]]
W=[[F(1),F(1)],[F(0),F(-1)]]
def field(Y):
 # F(Y)=X+Y_12 W. The field is nonlinear after multiplication by Y.
 out=[mscale(W,a[0][1]) for a in Y];out[0]=madd(out[0],X);return out

def matrix_check(U,V):
 p=5;truth=pzero(p);truth[0]=I
 for n in range(p):truth[n+1]=mscale(pmul(field(truth),truth)[n],F(1,n+1))
 results=[]
 for reverse in [False,True]:
  Y=pzero(p);Y[0]=I;Z=pzero(p)
  for u,v in zip(U,V):
   fy=field(Y);hfy=[zmat()]+fy[:-1];Z=padd(pscale(Z,u),hfy);arg=pscale(Z,v)
   E=pzero(p);E[0]=I;power=E.copy()
   for k in range(1,p+1):power=pmul(power,arg);E=padd(E,pscale(power,F(1,factorial(k))))
   Y=pmul(Y,E) if reverse else pmul(E,Y)
  residual=[madd(Y[n],mscale(truth[n],-1)) for n in range(p+1)]
  bad=[n for n,R in enumerate(residual) if R!=zmat()]
  if not reverse:require((bad and min(bad)==5), 'verify.py:135: check failed')
  else:require((bad and min(bad)<5), 'verify.py:136: check failed')
  results.append({'multiplication':'reversed negative control' if reverse else 'source convention','first_nonzero_local_error_degree':min(bad),
                  'leading_error_matrix':[[str(x) for x in row] for row in residual[min(bad)]]})
 return results

if __name__=='__main__':
 report={'status':'PASS','checks':identity_checks(),'emitted_polynomials_on_direct_series':certificate_example_checks()};U,V,p4=p4_checks();report['rational_method']=p4;report['noncommuting_matrix_check']=matrix_check(U,V)
 print(json.dumps(report,indent=2))
