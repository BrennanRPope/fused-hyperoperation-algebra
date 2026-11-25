import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S1,S2,S3,S4="#257de0","#f47230","#00b894","#6a4ce2"
INK,MUTE,GRID="#191919","#666666","#e6e6e6"
def setup():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.edgecolor":MUTE,"axes.labelcolor":INK,
        "xtick.color":MUTE,"ytick.color":MUTE,"axes.spines.top":False,"axes.spines.right":False,
        "axes.grid":True,"grid.color":GRID,"grid.linewidth":0.8,"lines.linewidth":2,"figure.dpi":110,
        "axes.titlesize":12,"axes.titleweight":"bold","legend.frameon":False})
def save(fig,path):
    fig.tight_layout(); fig.savefig(path,dpi=150,facecolor="white"); plt.close(fig)
