"""Inline SVG figures for the q page. Colours come from CSS variables on the page."""
import math, json
F={}
def knot_svg(fx, W=260, H=240, pad=18, n=540, mirror=False, label=''):
    pts=[fx(2*math.pi*i/n) for i in range(n+1)]
    if mirror: pts=[(x,y,-z) for x,y,z in pts]
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    sx=(W-2*pad)/(max(xs)-min(xs)); sy=(H-2*pad)/(max(ys)-min(ys)); s=min(sx,sy)
    cx=(max(xs)+min(xs))/2; cy=(max(ys)+min(ys))/2
    P=[(W/2+(x-cx)*s, H/2-(y-cy)*s, z) for x,y,z in pts]
    segs=[]
    step=6
    for i in range(0,n,step):
        seg=P[i:i+step+1]; z=sum(p[2] for p in seg)/len(seg)
        ext=[P[(i+j)%n] for j in range(-3,step+4)]
        d='M'+' L'.join(f'{x:.1f},{y:.1f}' for x,y,_ in seg)
        e='M'+' L'.join(f'{x:.1f},{y:.1f}' for x,y,_ in ext)
        segs.append((z,d,e))
    segs.sort()
    body=''.join(f'<path d="{d}" class="kb"/><path d="{e}" class="kn"/>' for z,d,e in segs)
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{label}">{body}</svg>'
tref=lambda t:(math.sin(t)+2*math.sin(2*t), math.cos(t)-2*math.cos(2*t), -math.sin(3*t))
fig8=lambda t:((2+math.cos(2*t))*math.cos(3*t),(2+math.cos(2*t))*math.sin(3*t),math.sin(4*t))
def hopf_svg(W=260,H=240):
    # two circles, linked: draw with gaps
    r=62; c1=(W/2-34,H/2); c2=(W/2+34,H/2)
    segs=[]
    n=240
    for k,(c,ph) in enumerate([(c1,0),(c2,math.pi)]):
        for i in range(0,n,4):
            a0=2*math.pi*i/n; a1=2*math.pi*(i+4)/n
            pt=lambda a:(c[0]+r*math.cos(a), c[1]+r*math.sin(a))
            xs=[pt(a0+(a1-a0)*j/4) for j in range(5)]
            ex=[pt(a0+(a1-a0)*j/4) for j in range(-2,7)]
            z=math.sin((a0+a1)/2)*(1 if k==0 else -1)
            d='M'+' L'.join(f'{x:.1f},{y:.1f}' for x,y in xs)
            e='M'+' L'.join(f'{x:.1f},{y:.1f}' for x,y in ex)
            segs.append((z,d,e))
    segs.sort()
    body=''.join(f'<path d="{d}" class="kb"/><path d="{e}" class="kn"/>' for z,d,e in segs)
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Hopf link">{body}</svg>'
F['trefoilL']=knot_svg(tref,label='left-handed trefoil')
F['trefoilR']=knot_svg(tref,mirror=True,label='right-handed trefoil')
F['fig8']=knot_svg(fig8,label='figure-eight knot')
F['hopf']=hopf_svg()

def cube_svg():
    # trefoil state cube: corners = states, layer by number of B-smoothings
    W,H=560,300
    loops={'000':3,'100':2,'010':2,'001':2,'110':1,'101':1,'011':1,'111':2}
    # layout by layer columns
    cols={0:['000'],1:['100','010','001'],2:['110','101','011'],3:['111']}
    pos={}
    for k,v in cols.items():
        for j,s in enumerate(v):
            x=70+k*140; y=H/2+(j-(len(v)-1)/2)*86
            pos[s]=(x,y)
    edges=''
    for s in pos:
        for i in range(3):
            if s[i]=='0':
                t=s[:i]+'1'+s[i+1:]
                (x1,y1),(x2,y2)=pos[s],pos[t]
                edges+=f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="ed"/>'
    nodes=''
    for s,(x,y) in pos.items():
        L=loops[s]
        nodes+=f'<circle cx="{x}" cy="{y}" r="25" class="nd"/>'
        # small loops drawn inside
        for j in range(L):
            nodes+=f'<circle cx="{x-(L-1)*6+j*12}" cy="{y-5}" r="5" class="lp"/>'
        nodes+=f'<text x="{x}" y="{y+15}" class="tl">{s}</text>'
    heads=''.join(f'<text x="{70+k*140}" y="22" class="th">{n} state{"s" if n>1 else ""} · A^{3-2*k}</text>' for k,n in zip(range(4),[1,3,3,1]))
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="state cube of the trefoil">{heads}{edges}{nodes}</svg>'
F['cube']=cube_svg()

def chord_svg():
    W,H=560,300; x0,x1,y0=120,440,120
    X=lambda s:x0+(x1-x0)*(s+0.25)/1.5
    Y=lambda c:y0-c*1300
    arc='M'+' L'.join(f'{X(s):.1f},{Y(s*(s-1)/4):.1f}' for s in [i/100*1.3-0.15 for i in range(101)])
    g=f'<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y0}" class="ax"/>'
    g+=f'<path d="{arc}" class="arc"/>'
    g+=f'<line x1="{X(0)}" y1="{y0}" x2="{X(1)}" y2="{y0}" class="chord"/>'
    for s,lab,cl,dx in [(0,'rung 0: addition','tll',-4),(1,'rung 1: mult-add','tlr',4)]:
        g+=f'<circle cx="{X(s)}" cy="{y0}" r="6" class="dot"/><text x="{X(s)+dx}" y="{y0+24}" class="{cl}">{lab}</text>'
    g+=f'<circle cx="{X(.5)}" cy="{Y(-1/16)}" r="5" class="dot2"/><text x="{X(.5)}" y="{Y(-1/16)+24}" class="tl">half rung: −1/16</text>'
    g+=f'<text x="{X(.5)}" y="{y0-12}" class="tlc">q dial: a + b + k·ab, no cubic term</text>'
    g+=f'<text x="{X(.5)}" y="{y0-80}" class="tl" style="fill:var(--pru)">rung flow B_s</text>'
    g+=f'<text x="{W/2}" y="{H-30}" class="ts">across: weight of ab (k on the dial, s on the flow)</text>'
    g+=f'<text x="{W/2}" y="{H-12}" class="ts">down: weight of a²b + ab², which is s(s−1)/4 on the flow</text>'
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="rung flow arc and q dial chord">{g}</svg>'
F['chord']=chord_svg()

def ferrers_svg():
    # Franklin: partition 7+6+5+3 : staircase 7,6,5 (length 3), smallest part 3 -> equal => move? choose 8+7+6+2
    parts=[8,7,6,2]; W,H=560,230; r=11; gap=26
    g=''
    def draw(parts,ox,hl_row,hl_diag,title):
        out=f'<text x="{ox}" y="22" class="thl">{title}</text>'
        for i,p in enumerate(parts):
            for j in range(p):
                cls='dt'
                if i==hl_row: cls='dr'
                elif hl_diag and j==p-1 and i<hl_diag: cls='dd'
                out+=f'<circle cx="{ox+j*gap+r}" cy="{50+i*gap}" r="{r-2}" class="{cls}"/>'
        return out
    g+=draw(parts,20,3,0,'8 + 7 + 6 + 2 (4 parts)')
    g+=draw([9,8,6],300,None,2,'9 + 8 + 6 (3 parts)')
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Franklin involution on Ferrers diagrams">{g}</svg>'
F['ferrers']=ferrers_svg()

def paths_svg():
    # all lattice paths in a 2x2 box, area under path
    import itertools
    W,H=560,150; cell=24
    paths=sorted(set(itertools.permutations('RRUU')))
    g=''
    for idx,p in enumerate(paths):
        ox=14+idx*92; oy=110
        # area above path inside box: count cells below path? use number of inversions (U before R)
        inv=sum(1 for i in range(4) for j in range(i+1,4) if p[i]=='U' and p[j]=='R')
        # shade cells under path (inversion count cells)
        x=y=0; pts=[(ox,oy)]
        for c in p:
            if c=='R': x+=1
            else: y+=1
            pts.append((ox+x*cell,oy-y*cell))
        poly=' '.join(f'{a},{b}' for a,b in pts)
        # area between path and top-left: region above-left of path
        area_poly=' '.join(f'{a},{b}' for a,b in pts)+f' {ox+2*cell},{oy}'
        g+=f'<rect x="{ox}" y="{oy-2*cell}" width="{2*cell}" height="{2*cell}" class="bx"/>'
        g+=f'<polygon points="{area_poly}" class="ar"/>'
        g+=f'<polyline points="{poly}" class="pt"/>'
        g+=f'<text x="{ox+cell}" y="{oy+24}" class="tl">q^{inv}</text>'
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="six lattice paths in a 2 by 2 box">{g}</svg>'
F['paths']=paths_svg()

def tate_svg():
    W,H=560,300; cx,cy=150,150; q=0.62
    g=''
    for k in range(4):
        r=120*q**k
        g+=f'<circle cx="{cx}" cy="{cy}" r="{r:.1f}" class="{"rg" if k in (0,1) else "rgf"}"/>'
    # the fundamental annulus shaded between r0 and q r0
    g=f'<path d="M{cx+120},{cy} A120,120 0 1,0 {cx-120},{cy} A120,120 0 1,0 {cx+120},{cy} M{cx+120*q},{cy} A{120*q},{120*q} 0 1,1 {cx-120*q},{cy} A{120*q},{120*q} 0 1,1 {cx+120*q},{cy}Z" class="ann" fill-rule="evenodd"/>'+g
    a=0.7
    for k in range(4):
        r=105*q**k
        g+=f'<circle cx="{cx+r*math.cos(a):.1f}" cy="{cy-r*math.sin(a):.1f}" r="5" class="dot"/>'
    g+=f'<text x="{cx+105*math.cos(a)+12:.1f}" y="{cy-105*math.sin(a)-6:.1f}" class="tsl">1 + a</text>'
    g+=f'<text x="{cx+65*math.cos(a)+40:.1f}" y="{cy-65*math.sin(a)+16:.1f}" class="tsl">q(1 + a)</text>'
    # torus on right
    tx,ty=420,150
    g+=f'<ellipse cx="{tx}" cy="{ty}" rx="110" ry="58" class="rg"/>'
    g+=f'<path d="M{tx-50},{ty-2} Q{tx},{ty+26} {tx+50},{ty-2}" class="rg"/>'
    g+=f'<path d="M{tx-36},{ty+8} Q{tx},{ty-12} {tx+36},{ty+8}" class="rg"/>'
    g+=f'<ellipse cx="{tx-70}" cy="{ty+8}" rx="16" ry="44" class="rgf"/>'
    g+=f'<text x="{tx}" y="{ty+92}" class="tl">glue the two edges of the ring: a torus</text>'
    g+=f'<text x="{cx}" y="{H-6}" class="tl">the ring between radius 1 and q</text>'
    g+=f'<path d="M265,150 L305,150" class="arrw"/>'
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Tate gluing: ring glued into a torus">{g}</svg>'
F['tate']=tate_svg()

def dial_svg():
    W,H=560,190; x0,x1,y=50,510,92
    X=lambda q:x0+(x1-x0)*(q+1.5)/4
    g=f'<line x1="{x0}" y1="{y}" x2="{x1}" y2="{y}" class="ax"/>'
    pts=[(-1,'−2','parity','1,0,1,0,…'),(0,'−1','failure union','1,1,1,…'),(1,'0','counting','1,2,3,…'),(2,'1','mult-add','1,3,7,15,…')]
    for q,k,name,seq in pts:
        g+=f'<circle cx="{X(q)}" cy="{y}" r="8" class="dot"/>'
        g+=f'<text x="{X(q)}" y="{y-40}" class="th">{name}</text>'
        g+=f'<text x="{X(q)}" y="{y-20}" class="ts">[n] = {seq}</text>'
        g+=f'<text x="{X(q)}" y="{y+30}" class="tl">q = {q}</text>'
        g+=f'<text x="{X(q)}" y="{y+50}" class="ts">k = {k}</text>'
    g+=f'<text x="{W/2}" y="{H-8}" class="ts">k = q − 1</text>'
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="the q dial">{g}</svg>'
F['dial']=dial_svg()

def ising_svg():
    W,H=280,260; cx,cy,R=140,128,90; n=12
    spins=[1,1,1,-1,-1,1,1,1,-1,-1,-1,1]
    g=f'<circle cx="{cx}" cy="{cy}" r="{R}" class="rgf"/>'
    for i,s in enumerate(spins):
        a=2*math.pi*i/n-math.pi/2; x=cx+R*math.cos(a); y=cy+R*math.sin(a)
        if spins[i]!=spins[(i+1)%n]:
            b=a+math.pi/n; g+=f'<line x1="{cx+(R-18)*math.cos(b):.1f}" y1="{cy+(R-18)*math.sin(b):.1f}" x2="{cx+(R+18)*math.cos(b):.1f}" y2="{cy+(R+18)*math.sin(b):.1f}" class="wall"/>'
        g+=f'<circle cx="{x:.1f}" cy="{y:.1f}" r="9" class="{"up" if s>0 else "dn"}"/>'
    g+=f'<text x="{cx}" y="{cy-4}" class="tl">4 walls</text><text x="{cx}" y="{cy+16}" class="ts">always even</text>'
    return f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Ising ring with domain walls">{g}</svg>'
F['ising']=ising_svg()
json.dump(F,open('page/figures.json','w'))
print({k:len(v) for k,v in F.items()})
