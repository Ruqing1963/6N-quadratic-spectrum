import csv, numpy as np, matplotlib, os
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
_DATA=os.path.normpath(os.path.join(_HERE,"..","data"))
_FIG=os.path.normpath(os.path.join(_HERE,"..","figures"))
os.makedirs(_DATA,exist_ok=True); os.makedirs(_FIG,exist_ok=True)
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.3,"figure.dpi":150,"savefig.bbox":"tight"})

def rd(n): return list(csv.DictReader(open(os.path.join(_DATA,n))))

# ============ FIGURE 1 : n^2+a ============
fig,ax=plt.subplots(1,2,figsize=(9.4,4.0))
sc=rd("n2a_Cscan.csv")
a=np.array([int(r["a"]) for r in sc]); C=np.array([float(r["C_product"]) for r in sc]); mask=[r["wing_mask"] for r in sc]
colmap={"right_only":"#1f77b4","left_only":"#d62728","both_1to2":"#2ca02c"}
cols=[colmap[m] for m in mask]
ax[0].axhline(1.0,color="0.6",lw=0.8,ls="--")
ax[0].bar(a,C,color=cols,width=0.8)
import matplotlib.patches as mp
ax[0].legend(handles=[mp.Patch(color=colmap[k],label=k.replace("_"," ")) for k in colmap],fontsize=6.5,loc="upper left")
ax[0].set_xlabel(r"shift $a$ in $n^2+a$"); ax[0].set_ylabel(r"richness $C(a)=\prod_q g_a(q)$")
ax[0].set_title(r"(A) the $a$-spectrum: mask (colour) $\times$ richness",fontsize=9)
# (B) measured vs predicted C for representative a
rp=rd("n2a_representative.csv")
Cp=np.array([float(r["C_product"]) for r in rp]); Ce=np.array([float(r["C_empirical"]) for r in rp])
lim=[0.4,2.1]
ax[1].plot(lim,lim,"k-",lw=0.7,alpha=0.6)
ax[1].scatter(Cp,Ce,s=26,color="#1f77b4",zorder=3)
for r in rp:
    ax[1].annotate("a=%s"%r["a"],(float(r["C_product"]),float(r["C_empirical"])),fontsize=6,xytext=(3,-2),textcoords="offset points")
ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel(r"predicted $C(a)=\prod_q g_a(q)$"); ax[1].set_ylabel(r"measured $Q_a(N)/\!\int dt/\log(t^2{+}a)$")
ax[1].set_title("(B) Bateman--Horn holds across the family",fontsize=9)
fig.suptitle(r"Mask rotation and richness in $n^2+a$ on the $6N$ skeleton",fontsize=10)
fig.savefig(os.path.join(_FIG,"p32_fig1.pdf")); print("fig1 done")

# ============ FIGURE 2 : Euler corridor & Heegner ladder ============
fig2,ax2=plt.subplots(1,2,figsize=(9.4,4.0))
co=rd("euler_corridor.csv")
q=np.array([int(r["q"]) for r in co]); o41=np.array([int(r["omega_A41"]) for r in co]); o9=np.array([int(r["omega_A9"]) for r in co])
ax2[0].vlines(q-0.6,0,o41,color="#1f77b4",lw=2.2,label=r"$n^2+n+41$ (disc $-163$)")
ax2[0].vlines(q+0.6,0,o9,color="#d62728",lw=1.4,alpha=0.8,label=r"$n^2+n+9$ (disc $-35$)")
ax2[0].axvspan(0,40,color="#1f77b4",alpha=0.06)
ax2[0].annotate("clean corridor:\n$\\omega(q)=0$ for all $q\\leq37$",(20,1.4),fontsize=7,ha="center",color="#1f77b4")
ax2[0].set_xlim(0,80); ax2[0].set_ylim(0,2.4); ax2[0].set_yticks([0,1,2])
ax2[0].set_xlabel(r"prime $q$"); ax2[0].set_ylabel(r"$\omega(q)$ = # roots of $f\equiv0$")
ax2[0].set_title("(A) the Heegner corridor of Euler's polynomial",fontsize=9)
ax2[0].legend(fontsize=7,loc="upper right")
# (B) Heegner ladder C vs A
fa=rd("euler_family.csv")
A=np.array([int(r["A"]) for r in fa]); Ce=np.array([float(r["C_empirical"]) for r in fa]); h1=np.array([int(r["class_number_1"]) for r in fa])
Al=A[h1==1]; Cl=Ce[h1==1]; order=np.argsort(Al)
ax2[1].plot(Al[order],Cl[order],"o-",ms=6,color="#1f77b4",label="class number 1 (Heegner)")
ax2[1].scatter(A[h1==0],Ce[h1==0],s=40,marker="x",color="#d62728",label="not class number 1")
for r in fa:
    ax2[1].annotate(r["A"],(int(r["A"]),float(r["C_empirical"])),fontsize=6.5,xytext=(4,2),textcoords="offset points")
ax2[1].set_xlabel(r"$A$ in $n^2+n+A$"); ax2[1].set_ylabel(r"measured richness $C$")
ax2[1].set_title(r"(B) the class-number-1 ladder, up to $n^2+n+41$",fontsize=9)
ax2[1].legend(fontsize=7,loc="upper left")
fig2.savefig(os.path.join(_FIG,"p32_fig2.pdf")); print("fig2 done")
