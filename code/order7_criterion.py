"""Exact finite checks for the conditional order-seven criterion.

This checks enumeration, determining-system ranks, null directions, and the
three continuous target values. It does NOT establish automatic order-seven
transfer for Williamson coefficients.
"""
from fractions import Fraction as F
from itertools import permutations
from functools import lru_cache

@lru_cache(None)
def shuffle(a,b):
    if not a:return {b:1}
    if not b:return {a:1}
    out={}
    for w,n in shuffle(a[1:],b).items():out[a[:1]+w]=out.get(a[:1]+w,0)+n
    for w,n in shuffle(a,b[1:]).items():out[b[:1]+w]=out.get(b[:1]+w,0)+n
    return out
def require(condition, detail):
    if not condition:
        raise ValueError(detail)

def echelon(rows):
    piv={}
    for original in rows:
        row=list(original)
        for j,r in sorted(piv.items()):
            c=row[j]
            if c:row=[a-c*b for a,b in zip(row,r)]
        j=next((j for j,x in enumerate(row) if x),None)
        if j is not None:
            c=row[j];piv[j]=[x/c for x in row]
    return piv

def matrix(pattern):
    words=sorted(set(permutations(pattern)))
    index={w:i for i,w in enumerate(words)}
    rows=[]
    for x in sorted(set(pattern)):
        rows.append([F(w[0]==x) for w in words])
    for w in words:
        for j in range(1,len(w)):
            row=[F(0)]*len(words)
            for v,n in shuffle(w[:j],w[j:]).items():row[index[v]]=F(n)
            rows.append(row)
    piv=echelon(rows)
    free=[j for j in range(len(words)) if j not in piv]
    null=[]
    for j in free:
        v=[F(0)]*len(words);v[j]=F(1)
        for k,r in sorted(piv.items(),reverse=True):
            v[k]=-sum(r[l]*v[l] for l in range(k+1,len(words)))
        require(any(v) and all(sum(a*b for a,b in zip(row,v)) == 0 for row in rows),
                ("invalid null direction", pattern))
        null.append([str(x) for x in v])
    answer={'pattern':pattern,'words':[''.join(w) for w in words],
            'rank':len(piv),'columns':len(words),'null':null}
    selected={"aaabb":"aabab", "aabc":"abac"}.get(pattern)
    if selected:
        position=answer["words"].index(selected)
        coordinate=[F(i==position) for i in range(len(words))]
        augmented=len(echelon(rows+[coordinate]))
        require(augmented == len(piv)+1 == len(words),
                ("selected coordinate did not remove the defect", pattern))
        require(all(F(vector[position]) != 0 for vector in null),
                ("selected coordinate vanishes on defect", pattern))
        answer["selected_coordinate"]=selected
        answer["augmented_rank"]=augmented
    return answer

# Generate ordinary rooted trees independently as multisets of smaller trees.
tree_levels={1:[()]}
for n in range(2,8):
    smaller=[(w,t) for w,ts in tree_levels.items() for t in ts]
    found=[]
    def forest(start,left,chosen):
        if not left:
            found.append(tuple(chosen));return
        for j in range(start,len(smaller)):
            w,t=smaller[j]
            if w<=left:forest(j,left-w,chosen+[t])
    forest(0,n-1,[])
    tree_levels[n]=found
counts={w:len(ts) for w,ts in tree_levels.items()}
alphabet=[(weight,k) for weight,n in counts.items() for k in range(n)]
families=[]
def collect(start,left,chosen):
    if not left:
        families.append(tuple(chosen));return
    for j in range(start,len(alphabet)):
        w,_=alphabet[j]
        if w<=left:collect(j,left-w,chosen+[alphabet[j]])
collect(0,7,[])
patterns={}
for family in families:
    multiplicities=sorted([family.count(x) for x in set(family)],reverse=True)
    pattern=''.join(chr(97+j)*n for j,n in enumerate(multiplicities))
    patterns.setdefault(pattern,[]).append(family)
import json
require(counts == {1:1,2:1,3:2,4:4,5:9,6:20,7:48},
        ("ordinary-tree counts", counts))
require(len(families)==115, ("weight-seven family count",len(families)))
expected={"a":(48,1,1),"ab":(37,2,2),"aab":(13,3,3),"abc":(5,6,6),
          "aaab":(5,4,4),"aaaab":(2,5,5),"aaaaab":(1,6,6),
          "aaaaaaa":(1,1,1),"aaabb":(1,9,10),"aabc":(2,11,12)}
require(set(patterns)==set(expected), ("pattern coverage",set(patterns)))
records=[]
for pattern,family_list in sorted(patterns.items()):
    record=dict(matrix(pattern),family_count=len(family_list),examples=family_list[:2])
    require((record["family_count"],record["rank"],record["columns"])==expected[pattern],
            ("pattern rank/count",pattern,record))
    records.append(record)

def monomial_simplex(weights, multiplier=F(1)):
    total=0
    for weight in weights:
        total+=weight
        multiplier/=total
    return multiplier

targets=[monomial_simplex((1,1,2,1,2)),
         monomial_simplex((1,2,1,3)),
         monomial_simplex((1,2,1,3),F(1,2))]
require(targets==[F(1,280),F(1,84),F(1,168)],("continuous targets",targets))
six_words=sorted(set(permutations("aaaaab")))
six_rows=[[F(word[0]==letter) for word in six_words] for letter in "ab"]
for left,right in (("a","aaaab"),("aa","aaab"),("aaa","aab"),("aaaa","ab")):
    multiplicity=shuffle(tuple(left),tuple(right))
    six_rows.append([F(multiplicity.get(word,0)) for word in six_words])
require(len(echelon(six_rows))==6,"Selected aaaaab matrix is singular")
print(json.dumps({'status':'PASS','ordinary_tree_counts':counts,
 'family_count':len(families),'patterns':records,
 'selected_aaaaab_matrix':[[int(x) for x in row] for row in six_rows],
 'selected_aaaaab_rank':6,
 'continuous_targets':list(map(str,targets)),
 'scope':'Conditional sufficient criterion only; automatic order-seven transfer remains open.'},
 indent=2))
