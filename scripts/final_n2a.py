import numpy as np, csv, os, time
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
from math import isqrt
t0=time.time()

Nmax=3_000_000
s=np.ones(Nmax+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(Nmax)+1):
    if s[p]: s[p*p::p]=False
oddp=[int(q) for q in np.nonzero(s)[0] if q>2]
Pscan=[q for q in oddp if q<=1_500_000]
print("primes ready %.1fs (%d odd; %d for C-scan)"%(time.time()-t0,len(oddp),len(Pscan)))

def tonelli(D,q):
    D%=q
    if D==0: return 0
    if q%4==3: return pow(D,(q+1)//4,q)
    Q=q-1; S=0
    while Q%2==0: Q//=2; S+=1
    z=2
    while pow(z,(q-1)//2,q)!=q-1: z+=1
    M=S;c=pow(z,Q,q);t=pow(D,Q,q);R=pow(D,(Q+1)//2,q)
    while t!=1:
        i=0;tt=t
        while tt!=1: tt=tt*tt%q;i+=1
        b=pow(c,1<<(M-i-1),q); M=i;c=b*b%q;t=t*c%q;R=R*b%q
    return R

def Cproduct(a,plist):
    C=1.0
    for q in plist:
        D=(-a)%q
        if D==0: om=1
        else: om=2 if pow(D,(q-1)//2,q)==1 else 0
        C*=(1-om/q)/(1-1/q)
    return C

def sieve_count(a):
    is_comp=np.zeros(Nmax+1,dtype=bool); is_comp[(a%2)::2]=True
    for q in oddp:
        D=(-a)%q
        if D==0: is_comp[q::q]=True
        else:
            if pow(D,(q-1)//2,q)==1:
                r=tonelli(D,q); is_comp[r::q]=True; is_comp[(q-r)::q]=True
    n0=2 if a>=0 else max(2,isqrt(1-a)+1)
    is_f=~is_comp; is_f[:n0]=False
    hi=isqrt(max(0,Nmax-a))+2
    for n in range(n0,min(hi,Nmax+1)):
        v=n*n+a
        if v<2: is_f[n]=False; continue
        pr=True
        for p in oddp:
            if p*p>v: break
            if v%p==0: pr=False; break
        if v%2==0 and v!=2: pr=False
        is_f[n]=pr
    nn=np.nonzero(is_f)[0]; Q=len(nn)
    vm6=((nn.astype(np.int64)**2+a)%6)
    right=int((vm6==1).sum()); left=int((vm6==5).sum())
    chan=sorted({int(x) for x in np.unique(nn%6)})
    tt=np.arange(n0,Nmax+1,dtype=np.float64); v=tt*tt+a; v=v[v>=3]
    BH=float(np.sum(1.0/np.log(v)))
    return Q,right,left,chan,(Q/BH if BH>0 else 0)

# representative table
reps=[-2,1,2,3,4,5,6,7,11,13]
with open(os.path.join(_DATA,"n2a_representative.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["a","a_mod6","Q","right_6Np1","left_6Nm1","right_over_left","channels_n_mod6","C_product","C_empirical"])
    print("\n  a a%6     Q     right   left   r/l   channels      Cprod  Cemp")
    for a in reps:
        Q,right,left,chan,Cemp=sieve_count(a); Cp=Cproduct(a,oddp)
        rl=("%.4f"%(right/left)) if left>0 else "inf"
        w.writerow([a,a%6,Q,right,left,rl,"|".join(map(str,chan)),"%.4f"%Cp,"%.4f"%Cemp])
        print(" %3d  %d  %7d %7d %7d %6s  %-12s %.4f %.4f"%(a,a%6,Q,right,left,rl,str(chan),Cp,Cemp))
print("rep table done %.1fs"%(time.time()-t0))

# C(a) richness scan (product only, fast)
with open(os.path.join(_DATA,"n2a_Cscan.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["a","a_mod6","C_product","wing_mask"])
    mask={0:"right_only",1:"both_1to2",2:"left_only",3:"right_only",4:"both_1to2",5:"left_only"}
    for a in range(-3,49):
        if a==0: continue
        Cp=Cproduct(a,Pscan)
        w.writerow([a,a%6,"%.4f"%Cp,mask[a%6]])
print("C-scan done %.1fs"%(time.time()-t0))
print("TOTAL %.1fs"%(time.time()-t0))
