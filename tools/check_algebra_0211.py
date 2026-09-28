"""Independent NumPy checks; these are not tests of the installed R package."""
import json, sys
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag
rng=np.random.default_rng(20260922)
errors={'score_relative':0.,'covariance_relative':0.,'half_relative':0.,'inv_half_relative':0.,'c_update_relative':0.}
counts={'scores':0,'roots':0,'updates':0}
def err(a,b):
 return float(np.max(np.abs(a-b),initial=0)/max(1.,float(np.max(np.abs(b),initial=0))))
def score(A,Q,ridge):
 d,V=np.linalg.eigh((Q+Q.T)/2)
 scale=np.diag(Q).mean()
 if scale<=0: return -np.inf
 d=np.maximum(d+ridge*scale,1e-10)
 return np.sum(np.sum(V*(A@V),axis=0)/d)
for n,ps in [(1,[1,3,7]),(4,[9,2,5]),(40,[2,3,4]),(12,[20,4,25])]:
 X=[rng.normal(size=(n,p))+rng.normal(size=p) for p in ps]
 XX=np.concatenate(X,axis=1); S=XX.T@XX/n
 S0=block_diag(*[x.T@x/n for x in X])
 for r in [1,2,4]:
  for dependent in [False,True]:
   L=rng.normal(size=(sum(ps),r))
   if dependent and r>1:L[:,-1]=L[:,0]
   blocks=np.split(L,np.cumsum(ps)[:-1]); T=[x@b for x,b in zip(X,blocks)]
   A=sum(T).T@sum(T)/n; Q=sum(t.T@t for t in T)/n
   Ad=L.T@S@L; Qd=L.T@S0@L
   np.testing.assert_allclose(A,Ad,rtol=1e-12,atol=1e-10)
   np.testing.assert_allclose(Q,Qd,rtol=1e-12,atol=1e-10)
   for ridge in [1e-8,1e-4]:
    x,y=score(A,Q,ridge),score(Ad,Qd,ridge)
    e=abs(x-y)/max(1,abs(y));errors['score_relative']=max(errors['score_relative'],e)
    assert e<2e-6,(n,ps,r,dependent,e)
    counts['scores']+=1
for n,p in [(5,17),(13,50),(30,12),(1,10)]:
 for zero in [False,True]:
  X=rng.normal(size=(n,p));X-=X.mean(0)
  X[:,0]=0
  if zero:X[:]=0
  u,d,vt=np.linalg.svd(X,full_matrices=False);keep=d>0;Q=vt.T[:,keep];d=d[keep]**2/n
  S=X.T@X/n;rec=(Q*d)@Q.T
  errors['covariance_relative']=max(errors['covariance_relative'],err(rec,S))
  scale=np.diag(S).mean();scale=scale if scale>0 else 1
  # Positive ridges are well-conditioned; zero is tested algebraically below
  # against the same exact low-rank model to avoid eigensolver nullspace noise.
  for ridge in [1e-4,1e-2]:
   dd,V=np.linalg.eigh((S+S.T)/2);dd=np.maximum(dd+ridge*scale,1e-10)
   base=max(ridge*scale,1e-10); ds=np.maximum(d+ridge*scale,1e-10)
   for power,key in [(0.5,'half_relative'),(-0.5,'inv_half_relative')]:
    a=(Q*(ds**power-base**power))@Q.T+base**power*np.eye(p)
    b=(V*dd**power)@V.T;e=err(a,b);errors[key]=max(errors[key],e)
    assert e<2e-9,(n,p,zero,power,e)
    counts['roots']+=1
for pk,pl,a,b in [(8,7,3,4),(1,9,1,2),(7,6,7,6),(5,4,0,2)]:
 Qk=np.linalg.qr(rng.normal(size=(pk,max(1,a))))[0][:,:a]
 Ql=np.linalg.qr(rng.normal(size=(pl,max(1,b))))[0][:,:b]
 dk=rng.uniform(.1,2,size=a);dl=rng.uniform(.1,2,size=b)
 Sk=(Qk*dk)@Qk.T;Sl=(Ql*dl)@Ql.T
 # Include a cross-block component outside the retained bases; it must survive.
 S=rng.normal(size=(pk,pl));St=Qk.T@S@Ql;rem=S-Qk@St@Ql.T
 for shift in [.1,1.,4.]:
  T=rng.normal(size=(pk,pl));Pt=Qk.T@T@Ql
  Ct=(St+shift*Pt)/(dk[:,None]*dl+shift)
  C=T+rem/shift+Qk@(Ct-Pt)@Ql.T
  B=S+shift*T
  lhs=Sk@C@Sl+shift*C;e=err(lhs,B)
  errors['c_update_relative']=max(errors['c_update_relative'],e)
  assert e<1e-11
  counts['updates']+=1
result={'status':'passed','scope':'Independent NumPy algebra only; not R execution or end-to-end benchmarking','cases':counts,'max_errors':errors}
Path(sys.argv[1] if len(sys.argv) > 1 else 'algebra_results_0211.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
