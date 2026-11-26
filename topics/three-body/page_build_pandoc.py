"""DRAFT: same page as page_build.py, but converts the Markdown with pandoc (for machines without the markdown module)."""
import re, sys, subprocess, types
src, out = sys.argv[1], sys.argv[2]
code = open(src + "/page_build.py").read()
fake = types.ModuleType("markdown")
fake.markdown = lambda md, extensions=None: subprocess.run(["pandoc", "-f", "gfm", "-t", "html"], input=md, capture_output=True, text=True, check=True).stdout
sys.modules["markdown"] = fake
exec(compile(code, "page_build.py", "exec"), {"__name__": "__main__"})
