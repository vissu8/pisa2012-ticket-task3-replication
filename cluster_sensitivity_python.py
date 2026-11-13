"""Sensitivity of the cluster solution to the number of clusters K (Section 4.3), standard library only.
Input files (same folder): data2_cluster_kmedoids_k8.csv, SeqN_CP2.csv
The reported solution corresponds to K-medoids (PAM) on the five standardized clustering variables with Euclidean distance;
the K = 8 run below reproduces the reported cluster labels exactly (adjusted Rand index = 1).
Run: python3 cluster_sensitivity_python.py   (takes a few minutes)
"""
import csv,collections,math
from math import comb
D=list(csv.DictReader(open('data2_cluster_kmedoids_k8.csv',encoding='utf-8-sig')))
keys=('SCIBO2','SCDBO','BOBO','BOBOBO','DA')
rows=[tuple(int(d[k]) for k in keys) for d in D]; lab=[int(d['clustering']) for d in D]
pats=sorted(set(rows)); idx={p:i for i,p in enumerate(pats)}; n=len(pats); w=[0]*n
for r in rows: w[idx[r]]+=1
mu=[sum(r[j] for r in rows)/len(rows) for j in range(5)]
sd=[math.sqrt(sum((r[j]-mu[j])**2 for r in rows)/(len(rows)-1)) for j in range(5)]
Z=[tuple((x-m)/s for x,m,s in zip(p,mu,sd)) for p in pats]
Dm=[[math.sqrt(sum((a-b)**2 for a,b in zip(Z[i],Z[j]))) for j in range(n)] for i in range(n)]
def pam(k):
    med=[min(range(n),key=lambda i:sum(w[j]*Dm[i][j] for j in range(n)))]
    while len(med)<k:
        near=[min(Dm[m][j] for m in med) for j in range(n)]
        med.append(max((c for c in range(n) if c not in med),key=lambda c:sum(w[j]*max(near[j]-Dm[c][j],0) for j in range(n))))
    cost=lambda M:sum(w[j]*min(Dm[m][j] for m in M) for j in range(n))
    cur=cost(med); improved=True
    while improved:
        improved=False
        for mi in range(k):
            for c in range(n):
                if c in med: continue
                M=med[:]; M[mi]=c; cc=cost(M)
                if cc<cur-1e-9: med,cur,improved=M,cc,True
    return [min(range(k),key=lambda t:Dm[med[t]][j]) for j in range(n)]
def ari(a,b):
    N=len(a); ct=collections.Counter(zip(a,b)); ra=collections.Counter(a); rb=collections.Counter(b)
    s=sum(comb(v,2) for v in ct.values()); sa=sum(comb(v,2) for v in ra.values()); sb=sum(comb(v,2) for v in rb.values()); e=sa*sb/comb(N,2)
    return (s-e)/((sa+sb)/2-e)
def v5(seq):
    t=seq.split()
    return (sum(x=='SCIBO2' for x in t),sum(x=='SCDBO' for x in t),sum(x.startswith('SCIBOBO') and not x.startswith('SCIBOBOBO') for x in t),sum(x.startswith('SCIBOBOBO') for x in t),sum(x in('SCD','SFI2','SCI2') for x in t))
S=[r for r in csv.DictReader(open('SeqN_CP2.csv',encoding='utf-8-sig')) if r['score']=='1']
cp=[(v5(r['seq']),float(r['CP2'])) for r in S]
def eta2(assign):
    g=collections.defaultdict(list)
    for p,y in cp: g[assign[idx[p]]].append(y)
    allv=[y for _,y in cp]; gm=sum(allv)/len(allv)
    return sum(len(v)*(sum(v)/len(v)-gm)**2 for v in g.values())/sum((y-gm)**2 for y in allv),g
res={}
print('distinct patterns: %d of %d correct responders; CP2 analyses n = %d'%(n,len(rows),len(cp)))
for k in (6,7,8,9,10):
    a=pam(k); res[k]=a; e,g=eta2(a); info=[]
    for c in range(k):
        mem=[i for i in range(n) if a[i]==c]; nn=sum(w[i] for i in mem)
        info.append((sum(w[i]*pats[i][4] for i in mem)/nn,nn,sum(g[c])/len(g[c]) if g.get(c) else float('nan')))
    info.sort(reverse=True)
    print('K=%d  eta2(CP2)=%.3f'%(k,e)+('  ARI vs reported K=8 labels = %.4f'%ari([a[idx[r]] for r in rows],lab) if k==8 else ''))
    print('   (mean N-CHECK, N, mean CP2): '+' | '.join('%.2f, %d, %.2f'%t for t in info))
print('ARI K=8 vs K=7: %.3f | K=8 vs K=9: %.3f'%(ari([res[8][idx[r]] for r in rows],[res[7][idx[r]] for r in rows]),ari([res[8][idx[r]] for r in rows],[res[9][idx[r]] for r in rows])))
