"""Embed results.json and the check count into the page template (index.html is self-contained)."""
import json, re, sys
tpl = open(sys.argv[1] if len(sys.argv) > 1 else "page_template.html").read()
data = json.load(open("results.json"))
keep = {k: data[k] for k in ["sim", "fano", "mu", "picture"]}
n = re.search(r"(\d+) PASS, 0 FAIL", open("checks_output.txt").read()).group(1)
open("index.html", "w").write(tpl.replace("__DATA__", json.dumps(keep, separators=(",", ":"))).replace("__CHECKS__", n))
print("wrote index.html with", n, "checks")
