"""Effect of the substitution on the five clustering variables and on the clusters (Online Resource 1, Table S6), standard library only.
Input files (same folder): data2_cluster_kmedoids_k8.csv, score1_raw_sequences.csv (raw action sequences of the 11,143 correct responders, in the order of score1_seqofseq.csv;
  one row per test taker, sequences separated by ' | ', a sequence without any selection is written as 'cancel'),
  score1_seqofseq.csv (SOS data after substitution), SeqN_CP2.csv, cancel_score1_notrank1.csv (56 key sequences and the similarity ranks of the keys).
Run: python3 raw_vs_substituted_python.py
"""
import csv,math,collections
raw_rows=[r['sequences'].split(' | ') for r in csv.DictReader(open('score1_raw_sequences.csv',encoding='utf-8-sig'))]
S=[x['sentence'].rstrip('.').split() for x in csv.DictReader(open('score1_seqofseq.csv',encoding='utf-8-sig'))]
keys=set(r['key'] for r in csv.DictReader(open('cancel_score1_notrank1.csv',encoding='utf-8-sig')))
assert len(raw_rows)==len(S)
total=sum(len(r) for r in raw_rows); empty=sum(sum(x=='cancel' for x in r) for r in raw_rows)
print('raw action sequences: %d | Cancel-only (no selection), removed: %d (%.1f%%)'%(total,empty,100*empty/total))
F=[[x for x in r if x!='cancel'] for r in raw_rows]
skip=[i for i,(f,s) in enumerate(zip(F,S)) if len(f)!=len(s)]
print('test takers whose number of sequences differs from the SOS (excluded):',len(skip))
use=[i for i in range(len(S)) if i not in set(skip)]
occ=exact=0
for i in use:
    for a,t in zip(F[i],S[i]):
        if t.endswith('BE'): continue
        occ+=1; exact+=(a in keys)
print('Cancel-concluding action sequences (n = %d test takers): %d | identical to a key: %d (%.1f%%) | substituted: %d (%.1f%%)'%(len(use),occ,exact,100*exact/occ,occ-exact,100*(occ-exact)/occ))
CH=('subway concession daily','subway full individual trip2','subway concession individual trip2')
raw=[sum(a in CH for a in f) for f in F]; fin=[sum(t in('SCD','SFI2','SCI2') for t in s) for s in S]
def corr(a,b):
    ma,mb=sum(a)/len(a),sum(b)/len(b); return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
print('N-CHECK without substitution vs after substitution: identical for %.1f%% | means %.3f vs %.3f | r = %.3f'%(100*sum(raw[i]==fin[i] for i in use)/len(use),sum(raw[i] for i in use)/len(use),sum(fin[i] for i in use)/len(use),corr([raw[i] for i in use],[fin[i] for i in use])))
Q=[r for r in csv.DictReader(open('SeqN_CP2.csv',encoding='utf-8-sig')) if r['score']=='1']
Sstr=[' '.join(s) for s in S]; Qstr=[r['seq'] for r in Q]
fw={};j=0
for i,s in enumerate(Sstr):
    if j<len(Qstr) and s==Qstr[j]: fw[i]=j;j+=1
bw={};j=len(Qstr)-1
for i in range(len(S)-1,-1,-1):
    if j>=0 and Sstr[i]==Qstr[j]: bw[i]=j;j-=1
ok=[i for i in fw if bw.get(i)==fw[i] and i not in set(skip)]
y=[float(Q[fw[i]]['CP2']) for i in ok]
print('CP2 linked through the order of the test takers (forward and backward matching agree): n = %d of %d'%(len(ok),len(Qstr)))
for lab,x in (('after substitution',fin),('without substitution',raw)):
    c=corr([x[i] for i in ok],y); print('N-CHECK %-22s r(CP2) = %.3f  R2 = %.3f'%(lab,c,c*c))


# ---- Five clustering variables computed from the raw sequences without the embedding, and the resulting clusters ----
# N-SCIBO2: Subway-Concession-Individual sequences with exactly one ticketing attempt (Buy-Ok clicks pressed in succession count once) and last trip choice Trip2
# N-SCDBO: Subway-Concession-Daily sequences with at least one Buy-Ok; N-BOBO / N-BOBOBO: Subway-Concession-Individual with two / three or more attempts
# N-CHECK: raw sequences identical to the three comparison actions
def attempts(q):
    t=q.split(); n=0
    for i,x in enumerate(t):
        if x=='buy' and i+1<len(t) and t[i+1]=='ok' and not (i>=2 and t[i-2]=='buy' and t[i-1]=='ok'): n+=1
    return n
def last_trip(q):
    t=q.split(); j=max(i for i,x in enumerate(t) if x=='buy' and i+1<len(t) and t[i+1]=='ok')
    for x in reversed(t[:j]):
        if x in('trip1','trip2','other'): return x
def rawvars(f):
    a=b=c=d=e=0
    for q in f:
        na=attempts(q)
        if q.startswith('subway concession individual'):
            if na==1 and last_trip(q)=='trip2': a+=1
            if na==2: c+=1
            if na>=3: d+=1
        if q.startswith('subway concession daily') and na>=1: b+=1
        if q in CH: e+=1
    return (a,b,c,d,e)
D=list(csv.DictReader(open('data2_cluster_kmedoids_k8.csv',encoding='utf-8-sig')))
V=('SCIBO2','SCDBO','BOBO','BOBOBO','DA'); finv=[tuple(int(d[k]) for k in V) for d in D]; lab=[int(d['clustering']) for d in D]
R=[rawvars(f) for f in F]
print()
for j,nm in enumerate(['N-SCIBO2','N-SCDBO','N-BOBO','N-BOBOBO','N-CHECK']):
    print('%-9s identical without and with the embedding: %.1f%%'%(nm,100*sum(R[i][j]==finv[i][j] for i in use)/len(use)))
mu=[sum(r[j] for r in finv)/len(finv) for j in range(5)]; sdv=[math.sqrt(sum((r[j]-mu[j])**2 for r in finv)/(len(finv)-1)) for j in range(5)]
zz=lambda p:tuple((x-m)/s for x,m,s in zip(p,mu,sdv))
cl=collections.defaultdict(list)
for p,l in zip(finv,lab): cl[l].append(p)
med={}
for l,rs in cl.items():
    c=collections.Counter(rs); med[l]=min(c,key=lambda p:sum(v*math.dist(zz(p),zz(q)) for q,v in c.items()))
pred={i:min(med,key=lambda t:math.dist(zz(R[i]),zz(med[t]))) for i in use}
print('cluster computed from the raw-sequence variables (nearest final medoid) equals the final cluster: %.1f%%'%(100*sum(pred[i]==lab[i] for i in use)/len(use)))
def eta2(lbl):
    g=collections.defaultdict(list)
    for i in ok: g[lbl[i]].append(float(Q[fw[i]]['CP2']))
    allv=[v for g_ in g.values() for v in g_]; gm=sum(allv)/len(allv)
    return sum(len(v)*(sum(v)/len(v)-gm)**2 for v in g.values())/sum((v-gm)**2 for v in allv)
print('CP2 variance associated with the clusters (n = %d): final %.3f | from raw-sequence variables %.3f'%(len(ok),eta2(lab),eta2(pred)))


# ---- Embedding alone: each Cancel-concluding sequence assigned to its first-ranked key, without the rules ----
import re
nr=list(csv.DictReader(open('cancel_score1_notrank1.csv',encoding='utf-8-sig')))
perseq=collections.defaultdict(list)
for r_ in nr: perseq[r_['seq']].append(r_)
keylist=sorted(keys); top1={}
for sq,rs in perseq.items():
    present={r_['key'] for r_ in rs}
    top1[sq]=[k for k in keylist if k not in present][0] if int(rs[0]['allocation'])==1 else [r_['key'] for r_ in rs if r_['rank']=='1'][0]
seen=collections.defaultdict(collections.Counter)
for i in use:
    for a,t in zip(F[i],S[i]): seen[a][t]+=1
tn={'trip1':'1','trip2':'2','other':'3+'}
def keycode(k):
    if k in seen: return seen[k].most_common(1)[0][0]
    n_=k.count('buy ok'); t_=re.findall(r'(trip1|trip2|other) buy ok',k); return ('SCIBOBO' if n_==2 else 'SCIBOBOBO')+tn[t_[-1]]
kc={k:keycode(k) for k in keylist}
tot=agree=sub=agree_sub=sub_c=agree_c=0; B=[]
for i in range(len(S)):
    if i in set(skip): B.append(finv[i]); continue
    toks=[]
    for a,t in zip(F[i],S[i]):
        if t.endswith('BE'): toks.append(t); continue
        c_=kc[top1[a]]; toks.append(c_); tot+=1; agree+=(c_==t)
        if a not in keys:
            sub+=1; agree_sub+=(c_==t)
            if t in ('SCD','SFI2','SCI2'): sub_c+=1; agree_c+=(c_==t)
    B.append((sum(x=='SCIBO2' for x in toks),sum(x=='SCDBO' for x in toks),sum(x.startswith('SCIBOBO') and not x.startswith('SCIBOBOBO') for x in toks),sum(x.startswith('SCIBOBOBO') for x in toks),sum(x in('SCD','SFI2','SCI2') for x in toks)))
print()
print('Embedding alone (first-ranked key, no rules): same code as the final assignment for %.1f%% of all Cancel-concluding sequences, for %.1f%% of the %d substituted sequences, and for %.1f%% of the %d substituted sequences that the final solution codes as comparison actions'%(100*agree/tot,100*agree_sub/sub,sub,100*agree_c/sub_c,sub_c))
for j,nm_ in enumerate(['N-SCIBO2','N-SCDBO','N-BOBO','N-BOBOBO','N-CHECK']):
    print('%-9s identical, embedding alone vs final: %.1f%%'%(nm_,100*sum(B[i][j]==finv[i][j] for i in use)/len(use)))
pb={i:min(med,key=lambda t:math.dist(zz(B[i]),zz(med[t]))) for i in use}
print('cluster from the embedding-alone variables equals the final cluster: %.1f%% (rules alone without the embedding: %.1f%%)'%(100*sum(pb[i]==lab[i] for i in use)/len(use),100*sum(pred[i]==lab[i] for i in use)/len(use)))
