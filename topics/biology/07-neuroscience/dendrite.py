"""Steady-state passive dendrite: soma + apical trunk + one oblique branch, compartmental."""
import numpy as np
def build(n_trunk=60,n_obl=25,obl_at=20,g_ax=30.0,g_leak=0.02,g_soma=0.3):
    n=1+n_trunk+n_obl; G=np.zeros((n,n))
    def link(a,b,g): G[a,a]+=g; G[b,b]+=g; G[a,b]-=g; G[b,a]-=g
    G[0,0]+=g_soma
    for m in range(1,n_trunk+1):
        G[m,m]+=g_leak; link(m-1,m,g_ax)
    for m in range(n_trunk+1,n):
        G[m,m]+=g_leak; link(obl_at if m==n_trunk+1 else m-1,m,g_ax*0.5)
    return G
def respond(G,inputs,s=0):
    """inputs: list of (node, conductance, reversal relative to rest). Returns somatic voltage."""
    A=G.copy(); b=np.zeros(len(G))
    for node,g,E in inputs: A[node,node]+=g; b[node]+=g*E
    return np.linalg.solve(A,b)[s]
def lens_constants(G,i,j,Ei,Ej,s=0):
    Z=np.linalg.inv(G)
    k=-Z[i,j]/(Ej*Z[s,i])-Z[i,j]/(Ei*Z[s,j]); q=Z[i,j]**2/(Ei*Ej*Z[s,i]*Z[s,j])
    return k,q
