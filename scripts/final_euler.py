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
print("primes ready %.1fs"%(time.time()-t0))

def tonelli(D,q):
    D%=q
    if D==0: return 0
    if q%4==3: return pow(D,(q+1)//4,q)
    Q=q-1;S=0
    while Q%2==0: Q//=2;S+=1
    z=2
    while pow(z,(q-1)//2,q)!=q-1: z+=1
    M=S;c=pow(z,Q,q);t=pow(D,Q,q);R=pow(D,(Q+1)//2,q)
    while t!=1:
        i=0;tt=t
        while tt!=1: tt=tt*tt%q;i+=1
        b=pow(c,1<<(M-i-1),q);M=i;c=b*b%q;t=t*c%q;R=R*b%q
    return R
inv2=lambda q:(q+1)//2

def omega(q,A):
    DD=(1-4*A)%q
    if DD==0: return 1
    return 2 if pow(DD,(q-1)//2,q)==1 else 0

def analyze(A):                       # A odd
    disc=1-4*A
    is_comp=np.zeros(Nmax+1,dtype=bool)
    first_div=None; Cp=1.0
    for q in oddp:
        DD=disc%q
        if DD==0:
            n0=(q-1)//2; is_comp[n0::q]=True; om=1
        else:
            if pow(DD,(q-1)//2,q)==1:
                m=tonelli(DD,q)
                for mm in (m,q-m):
                    is_comp[((mm-1)*inv2(q))%q::q]=True
                om=2
            else: om=0
        if om>0 and first_div is None: first_div=q
        Cp*=(1-om/q)/(1-1/q)
    is_f=~is_comp; is_f[0]=False
    for n in range(1,isqrt(Nmax)+2):
        v=n*n+n+A; pr=True
        for p in oddp:
            if p*p>v: break
            if v%p==0: pr=False; break
        if v%2==0 and v!=2: pr=False
        if n<len(is_f): is_f[n]=pr
    nn=np.nonzero(is_f)[0]; Q=len(nn)
    vm6=((nn.astype(np.int64)**2+nn+A)%6)
    right=int((vm6==1).sum()); left=int((vm6==5).sum())
    tt=np.arange(1,Nmax+1,dtype=np.float64); v=tt*tt+tt+A
    BH=float(np.sum(1.0/np.log(v)))
    return disc,first_div,Q,right,left,Cp,Q/BH

# family: class-number-1 (Heegner, odd A) and non-lucky odd contrasts
h1={1,3,5,11,17,41}    # 4A-1 in {3,11,19,43,67,163}, class number 1
fam=[1,3,5,11,17,41,7,9,21,33]
with open(os.path.join(_DATA,"euler_family.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["A","disc_1_minus_4A","class_number_1","first_prime_divisor","Q","right_6Np1","left_6Nm1","C_product","C_empirical"])
    print("\n  A   disc  h=1  firstdiv     Q     right   left   Cprod  Cemp")
    for A in fam:
        disc,fd,Q,right,left,Cp,Cemp=analyze(A)
        w.writerow([A,disc,int(A in h1),fd,Q,right,left,"%.4f"%Cp,"%.4f"%Cemp])
        print(" %3d %5d   %d   %4s   %8d %7d %7d  %.3f %.3f"%(A,disc,A in h1,str(fd),Q,right,left,Cp,Cemp))
print("family done %.1fs"%(time.time()-t0))

# corridor: omega(q) for n^2+n+41 (lucky) and n^2+n+9 (non-lucky) over small q
with open(os.path.join(_DATA,"euler_corridor.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["q","omega_A41","omega_A9"])
    for q in [p for p in oddp if p<=150]:
        w.writerow([q,omega(q,41),omega(q,9)])

# summary
with open(os.path.join(_DATA,"quad_summary.csv"),"w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["parameter","value"])
    d41,fd41,Q41,r41,l41,Cp41,Ce41=analyze(41)
    for k,v in [("Nmax",Nmax),("euler_A41_disc",-163),("euler_A41_first_divisor",fd41),
        ("euler_A41_Q",Q41),("euler_A41_C_empirical","%.4f"%Ce41),
        ("euler_A41_right_over_left","%.4f"%(r41/l41)),
        ("corridor_all_q_le_37_omega","0 (n^2+n+41 dodges every prime below 41)")]:
        w.writerow([k,v])
print("TOTAL %.1fs"%(time.time()-t0))
