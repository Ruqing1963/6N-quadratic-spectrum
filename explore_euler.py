import numpy as np, time
from math import isqrt
t0=time.time()
Nmax=2_000_000
s=np.ones(Nmax+1,dtype=bool); s[:2]=False
for p in range(2,isqrt(Nmax)+1):
    if s[p]: s[p*p::p]=False
oddp=[int(q) for q in np.nonzero(s)[0] if q>2]
print("primes ready %.1fs"%(time.time()-t0))

def tonelli(D,q):
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

def inv2(q): return (q+1)//2          # inverse of 2 mod odd q

def analyze(A):
    disc=1-4*A                         # (2n+1)^2 = disc (mod q) for f=0
    is_comp=np.zeros(Nmax+1,dtype=bool)
    # q=2: f=n^2+n+A ; n^2+n even -> f = A mod2 ; if A even, all n give even f
    if A%2==0: is_comp[:]=True         # A even -> f even always (degenerate)
    Cprod=1.0; omega={}; first_div=None
    for q in oddp:
        DD=disc%q
        if DD==0:                       # q | disc -> double root m=0 -> 2n+1=0 -> n=(q-1)/2
            n0=(q-1)//2; is_comp[n0::q]=True; om=1
        else:
            leg=pow(DD,(q-1)//2,q)
            if leg==1:
                m=tonelli(DD,q); 
                for mm in (m,q-m):
                    n0=((mm-1)*inv2(q))%q; is_comp[n0::q]=True
                om=2
            else: om=0
        if q<=60: omega[q]=om
        if om>0 and first_div is None: first_div=q
        Cprod*=(1-om/q)/(1-1/q)
    is_f=~is_comp; is_f[:1]=False
    for n in range(1,isqrt(Nmax)+2):
        v=n*n+n+A; pr=True
        for p in oddp:
            if p*p>v: break
            if v%p==0: pr=False; break
        if v%2==0 and v!=2: pr=False
        if n< len(is_f): is_f[n]=pr
    nn=np.nonzero(is_f)[0]; Q=len(nn)
    vm6=((nn.astype(np.int64)**2+nn+A)%6)
    right=int((vm6==1).sum()); left=int((vm6==5).sum()); other=int(((vm6!=1)&(vm6!=5)).sum())
    tt=np.arange(1,Nmax+1,dtype=np.float64); v=tt*tt+tt+A
    BH=float(np.sum(1.0/np.log(v))); Cemp=Q/BH
    return A,disc,Q,right,left,other,Cprod,Cemp,first_div,omega

print("\n n^2+n+A family (Euler lucky numbers have class number 1):")
print(" A    disc    Q       right   left  r/l   first q|f   C_prod   C_emp")
for A in [41,17,11,5,3,2,40,6,1]:
    A_,disc,Q,right,left,other,C,Cemp,fd,om=analyze(A)
    rl=("%.3f"%(right/left)) if left>0 else "inf"
    print(" %2d  %5d  %8d  %7d %7d %5s   %4s     %7.4f  %7.4f%s"%(
        A_,disc,Q,right,left,rl,str(fd),C,Cemp," other=%d"%other if other else ""))

print("\nomega(q) for n^2+n+41 (the corridor: 0 = prime never divisible by q):")
A_,disc,Q,right,left,other,C,Cemp,fd,om=analyze(41)
print("  ".join("q=%d:%d"%(q,om[q]) for q in sorted(om)))
print("first prime that can divide n^2+n+41:",fd)
print("\nTOTAL %.1fs"%(time.time()-t0))
