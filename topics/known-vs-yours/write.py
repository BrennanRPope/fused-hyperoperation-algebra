import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rows = json.load(open("rows.json")); OUT="/mnt/project-files/math-theory/known-vs-yours"
import os; os.makedirs(OUT, exist_ok=True)
G=["Known","Yours","Unsure"]; D=["2026-11-27","2026-11-28"]
c={(g,d):sum(1 for r in rows if r[2]==g and r[5]==d) for g in G for d in D}
fr=sum(1 for r in rows if r[2]=="Known" and r[3].startswith("your framing"))
fig,ax=plt.subplots(figsize=(8,3.4))
cols={"2026-11-27":"#8da0cb","2026-11-28":"#fc8d62"}
for i,g in enumerate(G):
    left=0
    for d in D:
        v=c[(g,d)]; ax.barh(i,v,left=left,color=cols[d],edgecolor="white",label=("Nov 27 prior-art audit" if d==D[0] else "Added Nov 28") if i==0 else None)
        if v: ax.text(left+v/2,i,str(v),ha="center",va="center",color="white",fontweight="bold")
        left+=v
    ax.text(left+0.8,i,f"{left} total",va="center")
ax.set_yticks(range(3)); ax.set_yticklabels(["Known\n(%d are your framing\nof known math)"%fr,"Yours\n(new result or new link)","Unsure\n(needs a lit check)"])
ax.invert_yaxis(); ax.set_xlim(0,max(sum(c[(g,d)] for d in D) for g in G)+10); ax.set_xlabel("items")
ax.set_title("DRAFT · Known vs yours, running tally (Nov 28 2026)"); ax.legend(loc="lower right",frameon=False)
for s in ["top","right"]: ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+"/fig1-tally.png",dpi=150)
def tab(g):
    L=["| # | Item | What's known, or exactly what's yours | Where | Added |","|---|---|---|---|---|"]
    for r in rows:
        if r[2]==g: L.append(f"| {r[0]} | {r[1]} | {r[3]} | `{r[4]}` | {r[5]} |")
    return "\n".join(L)
tot={g:sum(c[(g,d)] for d in D) for g in G}
md=f"""# DRAFT · Known vs yours: a running tally for Affine Composition Theory

*DRAFT, started Nov 28 2026. Brennan's ask (Nov 28): "we'll just have to keep a running tally". Built from the Nov 27
prior-art audit (`../prior-art-audit/`, 73 claims, about 110 searches) plus the newer write-ups' own "known vs ours"
sections and a short search today on the hyperoperation ladder. Not a claim of priority: "Yours" means nothing was found
after searching, never that nobody has it.*

## The answer first

- **{tot['Known']} known, {tot['Yours']} yours, {tot['Unsure']} unsure** ({len(rows)} items so far).
- The core law, the two-dial law and most physics rounds are known math or physics. **{fr} of the known items are your own
  framing of known results**, which is real work but not a new result.
- What's yours is mostly **links between fields** and **formal (Lean) work**, plus a few results: the fused hyperoperation
  ladder with its uniqueness proof, the Lee–Yang rule on the roundness dial, and the three-slit cap on the dial.
- **One correction:** the extended law's "new dial c" from the helix round is the repo's existing d dial (the same
  two-dial family, Coleman–McGuinness 1991). The helix reading on it is still open (U1).

![Counts per group](fig1-tally.png)

*Figure 1. Items per group. Blue came from the Nov 27 audit, orange was added Nov 28.*

## How to add to it

Append one row to the right table: a new id (K, Y or U plus the next number), one short line, the folder, the date.
Move a row between tables when a literature check settles it, and say in its line what settled it. Rows from the audit
keep its ids (A1, B2, M108 and so on); their sources are in `../prior-art-audit/README.md`.

Rules for sorting:
- **Known:** a source states it. If the project adds a reading of known pieces, the line starts "your framing".
- **Yours:** a real search found nothing, or the audit found the pieces but not the link. Say exactly what's new.
- **Unsure:** the search was thin or not done. These are the next literature checks.

## Known ({tot['Known']})

{tab('Known')}

## Yours ({tot['Yours']})

{tab('Yours')}

## Unsure ({tot['Unsure']})

{tab('Unsure')}

## Not yet sorted

The 145 write-ups were not read one by one. Rows come from the audit and from pages that have their own "known vs
ours" section. Pages without one (for example the chemistry, biology and NS planes pages) are not in the tally yet.

## Sources for today's search

- Hyperoperation background: https://en.wikipedia.org/wiki/Hyperoperation ; Mathlib `Mathlib.Data.Nat.Hyperoperation`.
- Distributivity of hyperoperations: Venkatesh, AMS Feature Column, https://mathvoices.ams.org/featurecolumn/?p=1510
- Bennett's commutative ladder: Jaramillo, https://www.arxiv.org/pdf/2012.04093
"""
open(OUT+"/README.md","w").write(md)
open("/mnt/project-files/math-theory/KNOWN-VS-YOURS.md","w").write(md.replace("](fig1-tally.png)","](known-vs-yours/fig1-tally.png)").replace("`../prior-art-audit/","`prior-art-audit/"))
print(tot, fr)
