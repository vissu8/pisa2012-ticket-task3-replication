"""Rank of the assigned key in the BERT (Max embedding) cosine-similarity ordering (Online Resource 1, Table S2).
Input: cancel_score1_notrank1.csv (columns: seq, key, score, rank, allocation, ...). One row per sequence-key pair;
'allocation' is the rank of the key finally assigned to the sequence after the domain-validity criteria, and that row is not in the file.
Run:  python3 substitution_audit_python.py   (file in the same folder)
"""
import csv,collections,math
rows=list(csv.DictReader(open('cancel_score1_notrank1.csv',encoding='utf-8-sig')))
alloc={}
for r in rows: alloc[r['seq']]=int(r['allocation'])
a=list(alloc.values()); n=len(a); m=sum(a)/n; sd=math.sqrt(sum((v-m)**2 for v in a)/n)
print('sequences: %d | mean rank of assigned key %.3f (population SD %.3f)'%(n,m,sd))
c=collections.Counter(a)
for lab,lo,hi in (('1',1,1),('2',2,2),('3',3,3),('4-5',4,5),('6-10',6,10),('>10',11,999)):
    k=sum(v for r_,v in c.items() if lo<=r_<=hi); print('rank %-5s n=%5d  %.1f%%'%(lab,k,100*k/n))
for t in (1,3,10): print('assigned key within top %d: %.1f%%'%(t,100*sum(v for r_,v in c.items() if r_<=t)/n))
keys=sorted(set(r['key'] for r in rows)); print('key sequences:',len(keys),'| with two or more Buy-Ok clicks:',sum(k.count('buy ok')>=2 for k in keys))


# ---- Assignment by token-level edit distance within the same admissibility criteria (Online Resource 1, Table S5) ----
# Admissible keys of a sequence = keys with check1 (Buy-Ok count) and check2 (final-ticket code) TRUE, plus the assigned key (its row is not in the file).
def lev(a,b):
    a=a.split(); b=b.split(); prev=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        cur=[i]
        for j,y in enumerate(b,1): cur.append(min(prev[j]+1,cur[j-1]+1,prev[j-1]+(x!=y)))
        prev=cur
    return prev[-1]
per=collections.defaultdict(list)
for r in rows: per[r['seq']].append(r)
cats=collections.Counter(); nadm=[]
for sq,rs in per.items():
    present={r['key']:r for r in rs}; ak=[k for k in keys if k not in present][0]
    adm=[k for k,r in present.items() if r['check1']=='TRUE' and r['check2']=='TRUE']+[ak]; nadm.append(len(adm))
    if len(adm)==1: cats['only one admissible key']+=1; continue
    d={k:lev(sq,k) for k in adm}; m=min(d.values()); best=[k for k in adm if d[k]==m]
    if len(best)==1: cats['edit distance unique, same key as BERT-based assignment' if best[0]==ak else 'edit distance unique, different key']+=1
    else: cats['edit distance tied, BERT-assigned key among nearest' if ak in best else 'edit distance tied, BERT-assigned key not among nearest']+=1
N=sum(cats.values()); print('\nadmissible keys per sequence: mean %.1f'%(sum(nadm)/N),sorted(collections.Counter(nadm).items()))
for k,v in cats.items(): print('%-62s %5d  %.1f%%'%(k,v,100*v/N))


# ---- Agreement with the similarity ordering compared with chance, for all sequences and for the sequences that are not themselves keys (RQ1) ----
# A key is always most similar to itself, so the 48 key sequences are excluded in the second line.
# Chance: if the k admissible keys of a sequence were random among the K keys, the expected rank of the first admissible key is (K+1)/(k+1).
# The number of admissible keys k counts the assigned key and the keys that satisfy the Buy-Ok-count and final-ticket criteria (check1, check2);
# the fare-class criterion could not be reconstructed from the file, so k is an upper bound and the chance level is conservative.
K=len(keys); kset=set(keys)
info={}
for sq,rs in per.items(): info[sq]=(int(rs[0]['allocation']),1+sum(1 for r in rs if r['check1']=='TRUE' and r['check2']=='TRUE'))
def comb_none(K,k,m):
    if K-k<m: return 0.0
    p_=1.0
    for i in range(m): p_*=(K-k-i)/(K-i)
    return p_
def agreement(seqs,label):
    n_=len(seqs); o=[info[q][0] for q in seqs]; kk=[info[q][1] for q in seqs]
    print('%s n=%d | mean rank %.2f (chance %.2f) | first-ranked key admissible %.1f%% (chance %.1f%%) | admissible key within the first three %.1f%% (chance %.1f%%)'%(label,n_,sum(o)/n_,sum((K+1)/(k+1) for k in kk)/n_,100*sum(v==1 for v in o)/n_,100*sum(k/K for k in kk)/n_,100*sum(v<=3 for v in o)/n_,100*sum(1-comb_none(K,k,3) for k in kk)/n_))
print()
agreement(list(info),'all sequences      ')
agreement([q for q in info if q not in kset],'non-key sequences  ')
