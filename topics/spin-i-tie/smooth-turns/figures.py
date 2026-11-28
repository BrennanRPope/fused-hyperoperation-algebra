# DRAFT figure for round 14
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,(a,b)=plt.subplots(1,2,figsize=(13,5.4))
t=np.linspace(0,2*np.pi,400); x=(1-np.exp(1j*t))/2
a.plot([0,1],[0,0],color='#9e9e9e',lw=4,label='real weights: the only road from 0 to 1')
a.plot(x.real,x.imag,color='#1565c0',lw=2.5,label='complex weights: lift = $e^{it}$, a smooth turn')
a.scatter([.5],[0],s=120,color='#c62828',zorder=5); a.annotate('x = 1/2: lift is 0\nthe state is erased,\nno way back',(.5,0),(.36,-.3),fontsize=10,color='#c62828',arrowprops=dict(arrowstyle='->',color='#c62828'))
for tt,lab in ((np.pi/2,'t = 90°'),(np.pi,'t = 180°: full flip, x = 1'),(3*np.pi/2,'t = 270°')):
    z=(1-np.exp(1j*tt))/2; a.scatter([z.real],[z.imag],color='#1565c0',zorder=5); a.annotate(lab,(z.real,z.imag),(z.real+.04,z.imag+.04),fontsize=9.5)
a.scatter([0],[0],color='k',zorder=5); a.annotate('x = 0: do nothing',(0,0),(-.22,.06),fontsize=9.5)
a.set_aspect('equal'); a.set_xlim(-.25,1.45); a.set_ylim(-.85,.65); a.set_xlabel('real part of x'); a.set_ylabel('imaginary part of x')
a.set_title('The XOR flip (k = -2) as a path of steps x',fontsize=11); a.legend(loc='lower right',fontsize=9)
T=np.linspace(0,4*np.pi,800)
b.plot(np.degrees(T),np.sin(T/2)**2,color='#1565c0',lw=2.5,label='flip odds $|x|^2 = \\sin^2(t/2)$')
b.axvline(360,color='#ef6c00',ls='--'); b.text(372,.42,"360°: the law's step is\nback to 'do nothing',\nbut a spin carries a\nminus sign here",fontsize=9.5,color='#ef6c00',bbox=dict(fc='white',ec='none'))
b.axvline(720,color='#2e7d32',ls='--'); b.text(728,.5,'720°:\nsign\ngone',fontsize=9.5,color='#2e7d32')
b.set_xlabel('turn angle t (degrees)'); b.set_ylabel('chance the bit has flipped'); b.set_xticks(range(0,721,90)); b.set_xlim(-20,790); b.set_ylim(-.02,1.12)
b.set_title('The odds follow the textbook qubit turn (half angle)',fontsize=11); b.legend(loc='upper left',fontsize=9.5)
fig.suptitle('DRAFT  round 14: going complex turns the law\'s flip into a smooth qubit turn about one axis; real weights cannot',fontsize=12,color='#b71c1c')
fig.tight_layout(); fig.savefig('st1_turns.png',dpi=130)
