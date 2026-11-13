"""Person-level analysis for Ticket Task 3 (standard library only).
Additional input for Table 9 and the checks of Tables 7-8: score1_seqofseq.csv (SOS data of the 11,143 correct responders, same row order as data2_cluster_kmedoids_k8.csv).
Inputs: SeqN_CP2.csv (score==1 -> 10,068 test takers) and data2_cluster_kmedoids_k8.csv (clustering variables and labels).
Cluster membership is reconstructed from the five-variable pattern, which maps to exactly one cluster in the clustering file.
Run:  python3 person_level_analysis_python.py
"""
import csv,glob,collections,math,re
U='./'  # folder containing the two CSV files
S=list(csv.DictReader(open(U+'SeqN_CP2.csv',encoding='utf-8-sig')))
D=list(csv.DictReader(open(U+'data2_cluster_kmedoids_k8.csv',encoding='utf-8-sig')))
corr=[r for r in S if r['score']=='1']
def vars_from_seq(s):
    t=s.split()
    return (sum(x=='SCIBO2' for x in t), sum(x=='SCDBO' for x in t), sum(x.startswith('SCIBOBO') and not x.startswith('SCIBOBOBO') for x in t),
            sum(x.startswith('SCIBOBOBO') for x in t), sum(x in ('SCD','SFI2','SCI2') for x in t))
pat2cl=collections.defaultdict(set); cnt2=collections.Counter()
for d in D:
    p=tuple(int(d[v]) for v in ('SCIBO2','SCDBO','BOBO','BOBOBO','DA')); pat2cl[p].add(int(d['clustering'])); cnt2[p]+=1

import math
lab={1:'I-1',2:'II-1',3:'III-1',4:'III-2',5:'II-3',6:'IV-1',7:'IV-2',8:'II-2'}
V5=('SCIBO2','SCDBO','BOBO','BOBOBO','DA')
rows=[]
for r in corr:
    p=vars_from_seq(r['seq']); c=next(iter(pat2cl[p]))
    rows.append(dict(cl=c,x=p,cp2=float(r['CP2']),seqn=float(r['seqn']),tt=float(r['totaltime']),fp=float(r['FirstPause']),lp=float(r['LastPause']),country=r['country']))
n=len(rows)
mean=lambda a:sum(a)/len(a)
def sd(a):
    m=mean(a); return math.sqrt(sum((x-m)**2 for x in a)/(len(a)-1))
print('N =',n)
paper={'I-1':(9.18,2.79,7.82,107.65),'II-1':(8.77,2.79,7.12,121.04),'II-2':(9.01,2.86,7.36,131.56),'II-3':(8.40,2.77,6.12,106.30),'III-1':(7.23,3.10,2.64,62.72),'III-2':(7.82,2.75,4.05,90.80),'IV-1':(7.49,3.01,4.07,102.07),'IV-2':(7.52,3.02,3.91,91.73)}
print('\ncluster | n | CP2 M(SD) | Seqn M | total time M | [paper: CP2 M(SD), Seqn, time]')
for c in sorted(lab,key=lambda k:['I-1','II-1','II-2','II-3','III-1','III-2','IV-1','IV-2'].index(lab[k])):
    g=[r for r in rows if r['cl']==c]
    print(f"{lab[c]:5s} n={len(g):5d} CP2 {mean([r['cp2'] for r in g]):.2f}({sd([r['cp2'] for r in g]):.2f}) Seqn {mean([r['seqn'] for r in g]):.2f} time {mean([r['tt'] for r in g]):.2f} | paper {paper[lab[c]]}")
print('overall CP2 %.2f(%.2f) Seqn %.2f time %.2f FP %.2f LP %.2f'%(mean([r['cp2'] for r in rows]),sd([r['cp2'] for r in rows]),mean([r['seqn'] for r in rows]),mean([r['tt'] for r in rows]),mean([r['fp'] for r in rows]),mean([r['lp'] for r in rows])))
# one-way ANOVA
gm=mean([r['cp2'] for r in rows]); groups={c:[r['cp2'] for r in rows if r['cl']==c] for c in lab}
ssb=sum(len(g)*(mean(g)-gm)**2 for g in groups.values()); ssw=sum(sum((x-mean(g))**2 for x in g) for g in groups.values())
F=(ssb/7)/(ssw/(n-8)); print('\nANOVA F(7,%d)=%.2f  eta2=%.4f   [paper: F=81.74, eta2=.054]'%(n-8,F,ssb/(ssb+ssw)))


import collections
n=len(rows)
lab={1:'I-1',2:'II-1',3:'III-1',4:'III-2',5:'II-3',6:'IV-1',7:'IV-2',8:'II-2'}
def inv(A):
    m=len(A); M=[row[:]+[1.0 if i==j else 0.0 for j in range(m)] for i,row in enumerate(A)]
    for c in range(m):
        p=max(range(c,m),key=lambda r:abs(M[r][c])); M[c],M[p]=M[p],M[c]; d=M[c][c]
        M[c]=[v/d for v in M[c]]
        for r in range(m):
            if r!=c and M[r][c]!=0:
                f=M[r][c]; M[r]=[a-f*b for a,b in zip(M[r],M[c])]
    return [row[m:] for row in M]
def ols(X,y,groups=None):
    p=len(X[0]); XtX=[[0.0]*p for _ in range(p)]; Xty=[0.0]*p
    for xi,yi in zip(X,y):
        for a in range(p):
            Xty[a]+=xi[a]*yi; xa=xi[a]; row=XtX[a]
            for b in range(a,p): row[b]+=xa*xi[b]
    for a in range(p):
        for b in range(a): XtX[a][b]=XtX[b][a]
    Ai=inv(XtX); beta=[sum(Ai[a][b]*Xty[b] for b in range(p)) for a in range(p)]
    res=[yi-sum(b*x for b,x in zip(beta,xi)) for xi,yi in zip(X,y)]
    rss=sum(e*e for e in res); df=len(y)-p; s2=rss/df
    se=[math.sqrt(s2*Ai[a][a]) for a in range(p)]
    out=dict(beta=beta,se=se,rss=rss,df=df,res=res,Ai=Ai)
    ybar=sum(y)/len(y); tss=sum((v-ybar)**2 for v in y); out['r2']=1-rss/tss
    if groups is not None:   # CR1 cluster-robust SE
        G=collections.defaultdict(lambda:[0.0]*p)
        for xi,e,g in zip(X,res,groups):
            for a in range(p): G[g][a]+=xi[a]*e
        meat=[[0.0]*p for _ in range(p)]
        for v in G.values():
            for a in range(p):
                for b in range(p): meat[a][b]+=v[a]*v[b]
        k=len(G); adj=(k/(k-1))*((len(y)-1)/(len(y)-p))
        # sandwich = Ai*meat*Ai
        T=[[sum(Ai[a][c]*meat[c][b] for c in range(p)) for b in range(p)] for a in range(p)]
        out['se_cr']=[math.sqrt(adj*sum(T[a][c]*Ai[c][a] for c in range(p))) for a in range(p)]
    return out
def pval(t): return math.erfc(abs(t)/math.sqrt(2))
y=[r['cp2'] for r in rows]; cty=[r['country'] for r in rows]
nchk=[r['x'][4] for r in rows]; lt=[math.log(max(r['tt'],0.1)) for r in rows]
so=[r['seqn']-r['x'][4] for r in rows]
print('min seqn_other',min(so),'| corr n_check vs cp2 =',end=' ')
def corr(a,b):
    ma,mb=sum(a)/len(a),sum(b)/len(b); return sum((x-ma)*(z-mb) for x,z in zip(a,b))/math.sqrt(sum((x-ma)**2 for x in a)*sum((z-mb)**2 for z in b))
print(round(corr(nchk,y),3))
print('correlations: nchk-logtime %.2f, nchk-seqn_other %.2f, logtime-seqn_other %.2f; cp2-logtime %.2f cp2-seqn_other %.2f'%(corr(nchk,lt),corr(nchk,so),corr(lt,so),corr(y,lt),corr(y,so)))
def fit(cols,names,groups=None):
    X=[[1.0]+[c[i] for c in cols] for i in range(n)]
    return ols(X,y,groups),['const']+names
def show(title,m,names,keys,robust=False):
    print('\n'+title+'  R2=%.4f'%m['r2'])
    for k in keys:
        i=names.index(k); b=m['beta'][i]; se=m['se_cr'][i] if robust else m['se'][i]
        sx=math.sqrt(sum((v-sum(c)/n)**2 for v in c)/(n-1)) if False else None
        print(f'  {k:10s} b={b:7.3f} SE={se:.3f} t={b/se:6.1f} p={pval(b/se):.1e}')
nm=['n_check'];m0,N0=fit([nchk],nm,cty)
m1,N1=fit([nchk,lt],['n_check','log_time'],cty)
m2,N2=fit([nchk,lt,so],['n_check','log_time','seqn_other'],cty)
show('M0: CP2 ~ N-CHECK',m0,N0,['n_check'],True); show('M1: + log(total time)',m1,N1,['n_check','log_time'],True); show('M2: + Seqn without N-CHECK',m2,N2,['n_check','log_time','seqn_other'],True)
# standardized beta and partial R2 for n_check in M2
sdy=math.sqrt(sum((v-sum(y)/n)**2 for v in y)/(n-1)); sdn=math.sqrt(sum((v-sum(nchk)/n)**2 for v in nchk)/(n-1))
print('  standardized beta (n_check, M0 / M2): %.3f / %.3f'%(m0['beta'][1]*sdn/sdy, m2['beta'][1]*sdn/sdy))
# partial R2 of n_check in M2: compare with model without n_check
mm,_=fit([lt,so],['log_time','seqn_other']); print('  M2 R2 %.4f vs without N-CHECK %.4f -> delta R2 = %.4f'%(m2['r2'],mm['r2'],m2['r2']-mm['r2']))
other=[[r['x'][j] for r in rows] for j in range(4)]
m3,N3=fit([nchk,lt,so]+other,['n_check','log_time','seqn_other','SCIBO2','SCDBO','BOBO','BOBOBO'],cty)
show('M3: + the other four clustering variables',m3,N3,['n_check','log_time','seqn_other','SCIBO2','SCDBO','BOBO','BOBOBO'],True)
# country fixed effects
cs=sorted(set(cty)); cd=[[1.0 if g==c else 0.0 for g in cty] for c in cs[1:]]
m4,N4=fit([nchk,lt,so]+cd,['n_check','log_time','seqn_other']+['c_'+c for c in cs[1:]],cty)
show('M4: M2 + country fixed effects',m4,N4,['n_check','log_time','seqn_other'],True)
# ---- cluster models
ref=3  # III-1 reference
dm={c:[1.0 if r['cl']==c else 0.0 for r in rows] for c in lab if c!=ref}
cn=['cl_'+lab[c] for c in dm]
c0,_=fit([dm[c] for c in dm],cn); c1,_=fit([dm[c] for c in dm]+[lt,so],cn+['log_time','seqn_other'])
base,_=fit([lt,so],['log_time','seqn_other'])
Fno=((sum(e*e for e in [v-sum(y)/n for v in y])-c0['rss'])/7)/(c0['rss']/c0['df'])
Fc=((base['rss']-c1['rss'])/7)/(c1['rss']/c1['df'])
print('\nCluster -> CP2 without controls: R2=%.4f (eta2)  F=%.2f'%(c0['r2'],Fno))
print('Cluster -> CP2 controlling log(total time) + Seqn w/o N-CHECK: R2 full=%.4f, controls only=%.4f; cluster partial eta2 = %.4f; F(7,%d)=%.2f'%(c1['r2'],base['r2'],(base['rss']-c1['rss'])/base['rss'],c1['df'],Fc))
# adjusted means at mean covariates
mlt=sum(lt)/n; mso=sum(so)/n
print('adjusted CP2 means (log_time, seqn_other at their means):')
names=['const']+cn+['log_time','seqn_other']
for c in sorted(lab,key=lambda k:['I-1','II-1','II-2','II-3','III-1','III-2','IV-1','IV-2'].index(lab[k])):
    b=c1['beta']; v=b[0]+(b[names.index('cl_'+lab[c])] if c!=ref else 0)+b[-2]*mlt+b[-1]*mso
    print(f'   {lab[c]:5s} {v:.2f}')
# does N-CHECK explain cluster effect? add N-CHECK to c1
c2,_=fit([dm[c] for c in dm]+[lt,so,nchk],cn+['log_time','seqn_other','n_check'])
print('\nCluster dummies + controls + N-CHECK: R2=%.4f ; N-CHECK b=%.3f'%(c2['r2'],c2['beta'][-1]))



# =====================================================================
# Position-based features (before / after the first cancellation message)
# Used for Table 9, the position-split model in Table 10, Cohen's d, and the R2 comparison in Section 4.5.
# The first cancellation message = first action sequence that contains a Buy-Ok click on a concession ticket
# (codes containing 'BO' but not 'BE').
# =====================================================================
CH=('SCD','SFI2','SCI2')
def is_trigger(t): return ('BO' in t) and not t.endswith('BE')
def next_cat(t):
    if t in ('SFI2BE','SFI1BE'): return 'Buy'
    if t in CH: return 'Compare'
    if t=='SCDBO': return 'Daily'
    if t.startswith('SCIBOBO'): return 'Repeat'
    if t.startswith('SCIBO'): return 'Same'
    return 'Other'
pos=[]
for r,row in zip([q for q in S if q['score']=='1'],rows):
    t=r['seq'].split(); i0=next(i for i,x in enumerate(t) if is_trigger(x))
    pos.append(dict(pre=sum(x in CH for x in t[:i0]),post=sum(x in CH for x in t[i0:]),first=(i0==0),nxt=next_cat(t[i0+1])))
print('\n--- Position of comparison actions (n = %d) ---'%n)
print('first purchase attempt in first action sequence: %.1f%%'%(100*sum(p['first'] for p in pos)/n))
print('compared before first message: %.1f%% ; after: %.1f%%'%(100*sum(p['pre']>0 for p in pos)/n,100*sum(p['post']>0 for p in pos)/n))
print('share of N-CHECK units before the first message: %.1f%%'%(100*sum(p['pre'] for p in pos)/sum(p['pre']+p['post'] for p in pos)))
order=['I-1','II-1','II-2','II-3','III-1','III-2','IV-1','IV-2']
cats=['Buy','Compare','Daily','Same','Repeat','Other']
print('\n(10,068 subsample) cluster | n | before M | after M | next action in percent: '+', '.join(cats))
for L in order+['Total']:
    idx=[i for i in range(n) if L=='Total' or lab[rows[i]['cl']]==L]
    g=len(idx)
    print('%-6s %5d  %.2f  %.2f  '%(L,g,sum(pos[i]['pre'] for i in idx)/g,sum(pos[i]['post'] for i in idx)/g)+'  '.join('%.1f'%(100*sum(pos[i]['nxt']==c for i in idx)/g) for c in cats))
pre=[p['pre'] for p in pos]; post=[p['post'] for p in pos]
m4,N4=fit([pre,post,lt,so],['pre','post','log_time','seqn_other'],cty)
show('M4: CP2 ~ N-CHECK before + after + log time + other action sequences',m4,N4,['pre','post','log_time','seqn_other'],True)
for lbl,cond in (('neither',lambda p:p['pre']==0 and p['post']==0),('after only',lambda p:p['pre']==0 and p['post']>0),('before only',lambda p:p['pre']>0 and p['post']==0),('both',lambda p:p['pre']>0 and p['post']>0)):
    g=[y[i] for i in range(n) if cond(pos[i])]
    print('  compared %-12s M=%.2f n=%d'%(lbl,mean(g),len(g)))
def cohen_d(a,b):
    na,nb=len(a),len(b); sp=math.sqrt(((na-1)*sd(a)**2+(nb-1)*sd(b)**2)/(na+nb-2)); return (mean(a)-mean(b))/sp
G=lambda names:[y[i] for i in range(n) if lab[rows[i]['cl']] in names]
print("\nCohen's d, I-1 vs III-1: %.2f ; I-1/II-1/II-2 vs III/IV: %.2f"%(cohen_d(G(['I-1']),G(['III-1'])),cohen_d(G(['I-1','II-1','II-2']),G(['III-1','III-2','IV-1','IV-2']))))
five=[[1.0]+list(r['x']) for r in rows]
print('R2 comparison: N-CHECK alone %.3f | five variables linear %.3f | eight clusters %.3f'%(m0['r2'],ols(five,y)['r2'],c0['r2']))


# =====================================================================
# Table 9 and checks of Tables 7-8 on all 11,143 correct responders
# (score1_seqofseq.csv is row-aligned with data2_cluster_kmedoids_k8.csv; row alignment is verified below)
# =====================================================================
SOS=[x['sentence'].rstrip('.').split() for x in csv.DictReader(open(U+'score1_seqofseq.csv',encoding='utf-8-sig'))]
assert len(SOS)==len(D)
aligned=sum(vars_from_seq(' '.join(t))==tuple(int(d[k]) for k in ('SCIBO2','SCDBO','BOBO','BOBOBO','DA')) for t,d in zip(SOS,D))
print('\n--- All correct responders (N = %d); SOS rows whose five variables match the clustering file: %d ---'%(len(SOS),aligned))
cl_all=[lab[int(d['clustering'])] for d in D]; N_all=len(SOS)
print('distinct SOS: %d | mean length %.2f (SD %.2f)'%(len(set(tuple(t) for t in SOS)),mean([len(t) for t in SOS]),sd([len(t) for t in SOS])))
pos_all=[]
for t in SOS:
    i0=next(i for i,x in enumerate(t) if is_trigger(x))
    pos_all.append(dict(pre=sum(x in CH for x in t[:i0]),post=sum(x in CH for x in t[i0:]),first=(i0==0),nxt=next_cat(t[i0+1])))
print('first purchase attempt in first action sequence: %.1f%%'%(100*sum(p['first'] for p in pos_all)/N_all))
print('compared before first message: %.1f%% ; after: %.1f%%'%(100*sum(p['pre']>0 for p in pos_all)/N_all,100*sum(p['post']>0 for p in pos_all)/N_all))
print('share of N-CHECK units before the first message: %.1f%%'%(100*sum(p['pre'] for p in pos_all)/sum(p['pre']+p['post'] for p in pos_all)))
print('\nTable 9: cluster | N | Seqn M (SD) | before M | after M | next action in percent: '+', '.join(cats))
for L in order+['Total']:
    idx=[i for i in range(N_all) if L=='Total' or cl_all[i]==L]; g=len(idx)
    sq=[len(SOS[i]) for i in idx]
    print('%-6s %5d %.2f (%.2f)  %.2f  %.2f  '%(L,g,mean(sq),sd(sq),sum(pos_all[i]['pre'] for i in idx)/g,sum(pos_all[i]['post'] for i in idx)/g)+'  '.join('%.1f'%(100*sum(pos_all[i]['nxt']==c for i in idx)/g) for c in cats))
g=[i for i in range(N_all) if cl_all[i]=='III-1']
print('III-1 dominant SOS (SCIBO2 SFI2BE): %.2f%%'%(100*sum(SOS[i]==['SCIBO2','SFI2BE'] for i in g)/len(g)))

# Sensitivity of the N-CHECK definition (Section 3.6): correlation with CP2 when the full-fare daily pass (SFD) is added to the comparison units
_S=[d for d in S if str(d['score']).strip() in ('1','1.0')]
_cp=[float(d['CP2']) for d in _S]
def _r(x,y):
    n=len(x); mx=sum(x)/n; my=sum(y)/n
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
for _name,_u in (('N-CHECK (SCD, SFI2, SCI2)',{'SCD','SFI2','SCI2'}),('N-CHECK + SFD',{'SCD','SFI2','SCI2','SFD'}),('SFD only',{'SFD'})):
    _x=[sum(t in _u for t in d['seq'].split()) for d in _S]
    print('%-26s r(CP2) = %.3f  (n = %d)'%(_name,_r(_x,_cp),len(_S)))
