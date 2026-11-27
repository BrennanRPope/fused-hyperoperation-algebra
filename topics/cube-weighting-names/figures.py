"""Figures for the cube-weighting naming page. Run: python3 figures.py"""
import itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BLUE, ORANGE, GREEN, PURPLE, GRAY, RED = "#2b6cb0", "#dd6b20", "#2f855a", "#6b46c1", "#4a5568", "#c53030"

# ---------- fig1: the cube, each merge costs one k ----------
fig, ax = plt.subplots(figsize=(8.2, 6.2))
P = {(0,0,0):(0,0),(1,0,0):(3,0),(0,1,0):(1.3,1.2),(0,0,1):(0,3),
     (1,1,0):(4.3,1.2),(1,0,1):(3,3),(0,1,1):(1.3,4.2),(1,1,1):(4.3,4.2)}
for a,b in itertools.combinations(P,2):
    if sum(abs(i-j) for i,j in zip(a,b))==1:
        ax.plot(*zip(P[a],P[b]),color="#a0aec0",lw=1.5,zorder=1)
lab = {(0,0,0):"1  (empty)",(1,0,0):"x",(0,1,0):"y",(0,0,1):"z",
       (1,1,0):"k·xy",(1,0,1):"k·xz",(0,1,1):"k·yz",(1,1,1):"k²·xyz"}
col = {0:GRAY,1:BLUE,2:ORANGE,3:RED}
for c,(px,py) in P.items():
    m=sum(c); ax.scatter(px,py,s=260,color=col[m],zorder=3)
    ax.text(px+0.18,py+0.18,lab[c],fontsize=13,color=col[m],zorder=4,fontweight="bold")
ax.text(5.7,3.9,"corner with m inputs\nweighs k^(m−1)",fontsize=12,color="black")
ax.text(5.7,2.6,"1 input:   weight 1",fontsize=11,color=BLUE)
ax.text(5.7,2.1,"2 inputs: weight k",fontsize=11,color=ORANGE)
ax.text(5.7,1.6,"3 inputs: weight k²",fontsize=11,color=RED)
ax.text(5.7,0.6,"x ∘ₖ y ∘ₖ z =\nsum of the 7 coloured corners",fontsize=11)
ax.set_title("Fig 1. The cube weighting: every merge of two corners costs one factor k",fontsize=12)
ax.set_xlim(-0.6,9.2); ax.set_ylim(-0.7,5.0); ax.axis("off")
fig.tight_layout(); fig.savefig("fig1-cube-weighting.png",dpi=150); plt.close(fig)

# ---------- fig2: one weighting, five names ----------
fig, ax = plt.subplots(figsize=(10,7))
def box(x,y,w,h,text,c,fs=10.5):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.02,rounding_size=0.08",
                 fc="white",ec=c,lw=2))
    ax.text(x,y,text,ha="center",va="center",fontsize=fs,color="black")
box(0,0.2,3.5,1.2,"corner weight k^(m−1)\n(our ∘ₖ folded over n inputs)",RED,12)
spokes = [
 (-4.2, 2.7,"Fuzzy measures\nSugeno λ-measure (1974)\nits Möbius / Harsanyi dividends\nare λ^(m−1)·∏gᵢ   (B2)",GREEN),
 ( 4.2, 2.7,"Formal groups\nmultiplicative formal group law\nx + y + βxy; the only polynomial\nones (A2, cube-weights)",BLUE),
 ( 4.6, 0.0,"K-theory of Grassmannians\nsingle-box stable Grothendieck\npolynomial G₁ (Fomin–Kirillov 1994,\nBuch 2002); β = k   (C1, C2)",PURPLE),
 (-4.6, 0.0,"Words and iterated sums\nweight-λ quasi-shuffle = mixable\nshuffle = Rota–Baxter weight λ\n(Hoffman 2000, Guo–Keigher 2000)  (D3, D4)",ORANGE),
 (-2.4,-3.0,"Multiple zeta values\nYamamoto's interpolation ζ^t (2013):\nmap Ψ with f = x/(1 − tx),\na block of m gets t^(m−1)   (E5)",RED),
 ( 2.4,-3.0,"Diehl–Ebrahimi-Fard–Tapia (2020)\ndeformed quasi-shuffle: weight k is\ntheir linear case f(t) = k t   (D2)",GRAY),
]
for x,y,t,c in spokes:
    box(x,y,4.3,1.5,t,c,9.4)
    ax.add_patch(FancyArrowPatch((0,0.2),(x*0.6,y*0.6+0.08),arrowstyle="-",color=c,lw=2,zorder=0))
ax.set_xlim(-6.8,6.8); ax.set_ylim(-4.1,3.8); ax.axis("off")
ax.set_title("Fig 2. The same weighting already has a name in six places (label in brackets = check)",fontsize=12)
fig.tight_layout(); fig.savefig("fig2-six-names.png",dpi=150); plt.close(fig)

# ---------- fig3: quasi-shuffle merges, ordered set partitions ----------
fig, ax = plt.subplots(figsize=(10,5.4))
def word(x0,y,blocks,c):
    for i,b in enumerate(blocks):
        h=0.42*len(b)
        ax.add_patch(FancyBboxPatch((x0+i*0.72,y),0.6,h,boxstyle="round,pad=0.01,rounding_size=0.05",fc="#edf2f7",ec=c,lw=1.6))
        for j,ch in enumerate(b):
            ax.text(x0+i*0.72+0.3,y+0.21+0.42*j,ch,ha="center",va="center",fontsize=12)
    return x0+len(blocks)*0.72
terms = [(["x"],["y"],["z"]),(["x"],["z"],["y"]),(["y"],["x"],["z"]),(["y"],["z"],["x"]),(["z"],["x"],["y"]),(["z"],["y"],["x"])]
x0=0.2
ax.text(0.2,4.75,"[x] ⋆ [y] ⋆ [z] at weight k  =",fontsize=13)
ax.text(0.2,4.2,"no merge (6 orders), weight 1:",fontsize=11,color=BLUE)
for i,t in enumerate(terms):
    word(0.2+i*2.45,3.3,t,BLUE)
ax.text(0.2,2.75,"one merge (6 words), weight k:",fontsize=11,color=ORANGE)
two=[(["x","y"],["z"]),(["z"],["x","y"]),(["x","z"],["y"]),(["y"],["x","z"]),(["y","z"],["x"]),(["x"],["y","z"])]
for i,t in enumerate(two):
    word(0.2+i*2.45,1.65,t,ORANGE)
ax.text(0.2,1.3,"everything merged (1 word), weight k²:",fontsize=11,color=RED)
word(0.4,-0.35,(["x","y","z"],),RED)
ax.text(1.6,0.3,"a stacked box is one letter made of m inputs; it costs k^(m−1),\nexactly the cube's corner weight (Fig 1). Keep only the one-box words and\nadd up 1 + [x] etc.: you get x ∘ₖ y ∘ₖ z  (check D4)",fontsize=10.5)
ax.set_xlim(0,15); ax.set_ylim(-0.6,5.1); ax.axis("off")
ax.set_title("Fig 3. Quasi-shuffle of weight k: ordered set partitions, each block of m weighs k^(m−1)",fontsize=12)
fig.tight_layout(); fig.savefig("fig3-quasi-shuffle-blocks.png",dpi=150); plt.close(fig)

# ---------- fig4: L_k does the same job twice ----------
fig, ax = plt.subplots(figsize=(10,5.2))
def node(x,y,t,c):
    ax.add_patch(FancyBboxPatch((x-1.7,y-0.45),3.4,0.9,boxstyle="round,pad=0.02,rounding_size=0.08",fc="white",ec=c,lw=2))
    ax.text(x,y,t,ha="center",va="center",fontsize=11)
def arrow(a,b,t,c,off=(0,0.25)):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=16,color=c,lw=1.8))
    ax.text((a[0]+b[0])/2+off[0],(a[1]+b[1])/2+off[1],t,ha="center",fontsize=10.5,color=c)
node(-3,1.6,"numbers with x ∘ₖ y",BLUE); node(3,1.6,"numbers with x + y",BLUE)
arrow((-1.25,1.6),(1.25,1.6),"L_k(x) = log(1 + kx)/k",BLUE)
node(-3,-1.0,"words with ⋆ of weight k",ORANGE); node(3,-1.0,"words with the shuffle",ORANGE)
arrow((-1.25,-1.0),(1.25,-1.0),"Ψ_{L_k}  (Hoffman's log)",ORANGE)
ax.text(0,0.45,"same series:  L_k(t) = t − k t²/2 + k² t³/3 − …   (checks E1, E4)",ha="center",fontsize=11)
ax.text(0,-2.2,"top row: one input at a time, order forgotten.  bottom row: inputs kept in time order.\nThe back arrows are (e^{kt} − 1)/k on both rows (E2, E3).",ha="center",fontsize=10)
ax.set_xlim(-5.2,5.2); ax.set_ylim(-2.7,2.4); ax.axis("off")
ax.set_title("Fig 4. Our flat coordinate L_k is the series in Hoffman's log map",fontsize=12)
fig.tight_layout(); fig.savefig("fig4-flat-coordinate-is-hoffman-log.png",dpi=150); plt.close(fig)
print("wrote 4 figures")
