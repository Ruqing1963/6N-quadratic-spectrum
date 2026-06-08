import numpy as np, time
from math import isqrt, log
t0=time.time()
Nmax=2_000_000

# primes up to Nmax
s=np.ones(Nmax+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(Nmax)+1):
    if s[p]: s[p*p::p]=False
primes=np.nonzero(s)[0]
oddp=[int(q) for q in primes if q>2]
print("primes<=Nmax=%d : %d  (%.1fs)"%(Nmax,len(primes),time.time()-t0))

def tonelli(D,q):                       # sqrt of D mod q (q odd prime, (D|q)=1)
    D%=q
    if D==0: return 0
    if q%4==3: return pow(D,(q+1)//4,q)
    Q=q-1; S=0
    while Q%2==0: Q//=2; S+=1
    z=2
    while pow(z,(q-1)//2,q)!=q-1: z+=1
    M=S; c=pow(z,Q,q); t=pow(D,Q,q); R=pow(D,(Q+1)//2,q)
    while t!=1:
        i=0; tt=t
        while tt!=1: tt=tt*tt%q; i+=1
        b=pow(c,1<<(M-i-1),q)
        M=i; c=b*b%q; t=t*c%q; R=R*b%q
    return R

def analyze(a):
    is_comp=np.zeros(Nmax+1,dtype=bool)
    is_comp[(a%2)::2]=True              # q=2: n^2+a even iff n=a (mod2)
    Cprod=1.0
    for q in oddp:
        D=(-a)%q
        if D==0:                        # q|a -> q|n^2+a iff q|n ; omega=1
            is_comp[q::q]=True; om=1
        else:
            leg=pow(D,(q-1)//2,q)
            if leg==1:
                r=tonelli(D,q); is_comp[r::q]=True; is_comp[(q-r)::q]=True; om=2
            else: om=0
        Cprod*=(1-om/q)/(1-1/q)
    # candidate n: n>=n0 so that n^2+a>=2
    n0=2
    if a<0: n0=max(2,isqrt(1-a)+1)
    is_f=~is_comp; is_f[:n0]=False
    # override small n (n^2+a <= Nmax) by direct primality
    lim=isqrt(Nmax) if Nmax-a>0 else 0
    for n in range(n0,isqrt(max(0,Nmax-a))+2):
        v=n*n+a
        if v<2: is_f[n]=False; continue
        pr=True
        for p in oddp:
            if p*p>v: break
            if v%p==0: pr=False; break
        if v%2==0 and v!=2: pr=False
        is_f[n]=pr
    nn=np.nonzero(is_f)[0]
    if len(nn)==0: return a,0,0,0,Cprod,[]
    val_mod6=( (nn.astype(np.int64)**2 + a) %6 )
    right=int((val_mod6==1).sum()); left=int((val_mod6==5).sum())
    other=int(((val_mod6!=1)&(val_mod6!=5)).sum())   # should be 0 for primes>3
    # which n mod6 channels are active
    chan=sorted(set(int(x) for x in np.unique(nn%6)))
    # Bateman-Horn integral
    tt=np.arange(n0,Nmax+1,dtype=np.float64); v=tt*tt+a
    BH=float(np.sum(1.0/np.log(v[v>=3])))
    Q=len(nn)
    return a,Q,right,left,Cprod,chan,other,Q/BH

print("\n  a  a%6   Q       right    left   r/l     channels(n mod6)  C_prod  C_emp(Q/int)")
for a in [1,7,2,3,4,5,6,-2,-1]:
    res=analyze(a)
    a_,Q,right,left,C,chan,other,Cemp=res
    rl = right/left if left>0 else float('inf')
    print(" %3d  %d  %7d  %7d %7d  %5.3f   %-16s  %.4f  %.4f%s"%(
        a_,a_%6,Q,right,left,rl,str(chan),C,Cemp," other=%d!"%other if other else ""))
print("\nTOTAL %.1fs"%(time.time()-t0))
