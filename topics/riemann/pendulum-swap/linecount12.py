"""linecount.py for n = 12, one frame and a slice of angles per process; writes one line per angle as it goes. DRAFT."""
import numpy as np, sys, json
sys.path.insert(0, "../rewind-dial"); sys.path.insert(0, "../knauf-field")
from knauf import chain
from oddpath import setup_odd
from oddrewind_dual import setup_frame
from frames import to_bonds
from linecount import counts
if __name__ == "__main__":
    frame, d0, d1, step = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); n = 12
    S = setup_odd(n) if frame == "spin" else setup_frame(to_bonds(np.log(chain(n).astype(float)), n), n)
    with open("lc12_%s_%d_%d.jsonl" % (frame, d0, d1), "a") as f:
        for deg in range(d0, d1, step):
            on, tot = counts(S, np.radians(deg)); f.write(json.dumps([deg, on, tot]) + "\n"); f.flush()
