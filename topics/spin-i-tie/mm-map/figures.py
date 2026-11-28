# DRAFT figure for round 13
import json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=json.load(open('results.json'))['two_bit_counts']
fig,(a,b)=plt.subplots(1,2,figsize=(13,5.2),gridspec_kw={'width_ratios':[1.25,1]})
rows=[('Frame: mixing is linear','touched','coin rule pins k = 0 (round 8)'),
      ('R1 Finiteness','no','the law says nothing'),
      ('R2 Local tomography','yes','free numbers compose by the law at k = +1'),
      ('R3 Equivalence of subspaces','no','the law says nothing'),
      ('R4 Symmetry (reversible turns)','next round','do the flips build smooth turns?'),
      ('R5 All measurements allowed','no','the law says nothing')]
col={'yes':'#2e7d32','touched':'#1565c0','no':'#9e9e9e','next round':'#ef6c00'}
a.axis('off'); a.set_xlim(0,1); a.set_ylim(0,len(rows)+0.8)
a.text(0,len(rows)+0.4,"Masanes & Mueller's rules vs the law",fontsize=13,weight='bold')
for i,(r,m,note) in enumerate(rows):
    y=len(rows)-i-0.4
    a.add_patch(plt.Rectangle((0,y-0.38),0.03,0.76,color=col[m]))
    a.text(0.05,y+0.12,r,fontsize=11,weight='bold'); a.text(0.05,y-0.25,note,fontsize=9.5,color='#444')
    a.text(0.98,y,m,fontsize=10.5,color=col[m],ha='right',weight='bold')
names=list(R); x=range(len(names))
b.bar([i-0.2 for i in x],[R[n]['N_AB actual'] for n in names],0.4,label='numbers two bits actually need',color='#455a64')
b.bar([i+0.2 for i in x],[R[n]['N_AB by law k=1'] for n in names],0.4,label='what the law at k = +1 predicts',color='#ffb300')
b.set_xticks(list(x)); b.set_xticklabels(names,fontsize=9.5)
b.set_ylabel('free numbers to pin down a two-bit state'); b.set_ylim(0,44); b.legend(loc='upper left',fontsize=9.5)
b.set_title('R2: only theories where the bars match are locally tomographic',fontsize=11)
fig.suptitle('DRAFT  round 13: local tomography (R2) is the law at k = +1; three rules lie outside it',fontsize=12,color='#b71c1c')
fig.tight_layout(); fig.savefig('mm1_map.png',dpi=130)
