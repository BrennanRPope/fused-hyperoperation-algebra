# DRAFT figure for round 16
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.linalg import expm
X=np.array([[0,1],[1,0]],complex); Z=np.diag([1,-1]).astype(complex); Y=-1j*Z@X
bloch=lambda psi: np.real([psi.conj()@P@psi for P in (X,Y,Z)])
fig=plt.figure(figsize=(13,5.6)); a=fig.add_subplot(1,2,1,projection='3d'); b=fig.add_subplot(1,2,2)
u,v=np.mgrid[0:2*np.pi:40j,0:np.pi:20j]; a.plot_wireframe(np.cos(u)*np.sin(v),np.sin(u)*np.sin(v),np.cos(v),color='#ddd',lw=.5)
psi0=np.array([1,1],complex)/np.sqrt(2); psi=psi0.copy(); legs=[(Z,1.4,'#1565c0','turn about Z'),(X,1.0,'#c62828','turn about X'),(Z,1.2,'#1565c0',None)]
pts=[bloch(psi)]
for P,ang,col,lab in legs:
    tr=[bloch(expm(-1j*s*P/2)@psi) for s in np.linspace(0,ang,60)]; tr=np.array(tr)
    a.plot(*tr.T,color=col,lw=3,label=lab); psi=expm(-1j*ang*P/2)@psi
a.scatter(*bloch(psi0),color='k',s=40); a.scatter(*bloch(psi),color='#2e7d32',s=60)
a.text(*(bloch(psi0)*1.2+np.array([0,0,.08])),'start',fontsize=9.5); e=bloch(psi); a.text(*(e*1.15+np.array([0,-.45,.12])),'reached',fontsize=9.5,color='#2e7d32')
for vec,lab in (((1.3,0,0),'X'),((0,1.3,0),'Y'),((0,0,1.3),'Z')): a.quiver(0,0,0,*vec,color='#555',arrow_length_ratio=.08); a.text(*(np.array(vec)*1.08),lab)
a.set_axis_off(); a.legend(loc='lower left',fontsize=9); a.set_title('Two flips, each turned smoothly by the law,\nreach any direction (here Z, then X, then Z,\nfrom a start on the equator)',fontsize=11)
b.axis('off'); b.set_xlim(0,1); b.set_ylim(0,6)
rows=[('X · X','+1','a flip (XOR)','#2e7d32'),('Z · Z','+1','a flip (XOR)','#2e7d32'),('X Z  vs  Z X','opposite signs','the one extra rule: they anticommute','#ef6c00'),
      ('(ZX) · (ZX)','-1','a quarter turn, not a flip','#c62828'),('(-i ZX) · (-i ZX)','+1','with i, the third axis Y is a flip again','#1565c0')]
b.text(0,5.6,'What each combination does when done twice',fontsize=12,weight='bold')
for i,(e_,val,note,col) in enumerate(rows):
    y=4.8-i*0.95; b.add_patch(plt.Rectangle((0,y-.35),.025,.7,color=col))
    b.text(.05,y+.1,e_,fontsize=12,family='monospace'); b.text(.55,y+.1,'= '+val,fontsize=12,color=col,weight='bold'); b.text(.05,y-.25,note,fontsize=9.5,color='#444')
fig.suptitle('DRAFT  round 16: two anticommuting XOR flips glue into full 3D spin, and the third axis needs i',fontsize=12,color='#b71c1c')
fig.tight_layout(); fig.savefig('s3_spin.png',dpi=130)
