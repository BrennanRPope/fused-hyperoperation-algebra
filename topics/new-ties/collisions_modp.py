"""Round 5: exact no-collision check for the round-3/4 example (walls ±3) to longer words, by arithmetic modulo a
large prime. If the maps are all different modulo P they are all different over the rationals, and the number of
different maps can never exceed the trace-monoid count; so 'count mod P = trace count' proves no collisions up to that
length (a collision would lower the count)."""
import io, contextlib, json
from fractions import Fraction as F
with contextlib.redirect_stdout(io.StringIO()):
    import planelaws as P
import union_ifs as U
PR=(1<<61)-1
def red(x): return x.numerator*pow(x.denominator,-1,PR)%PR
def canon(a,b,c,d):
    for t in (d,c,b,a):
        if t: inv=pow(t,-1,PR); return (a*inv%PR,b*inv%PR,c*inv%PR,d*inv%PR)
def run(q,n):
    G=[]
    for p,qq,ms in [(F(1),F(q),[F(1,5),F(3,20)]),(F(-1),F(-q),[F(1,4)])]:
        for s in ms:
            M=P.map_with(p,qq,s)[0]; G.append(tuple(red(x) for r in M for x in r))
    lv={(1,0,0,1)}; counts=[]
    for k in range(n):
        new=set()
        for (a,b,c,d) in lv:
            for (e,f,g,h) in G:
                new.add(canon((e*a+f*c)%PR,(e*b+f*d)%PR,(g*a+h*c)%PR,(g*b+h*d)%PR))
        lv=new; counts.append(len(lv))
    tc=U.trace_counts([(F(0),[F(0),F(1)]),(F(1),[F(0)])],n)
    return counts,tc
if __name__=="__main__":
    out={}
    for q in (2,3,12):
        c,t=run(q,14); out[q]=dict(counts=c,trace=t,match=c==t); print(q,c[-3:],t[-3:],c==t)
    json.dump(out,open("collisions_modp.json","w"),indent=1)
