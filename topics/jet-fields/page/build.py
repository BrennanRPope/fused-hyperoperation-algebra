"""Rebuild page/index.html: run figures.py first, then this."""
import os
d = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(d, 'template.html')).read()
t = t.replace('%%CUBE%%', open(os.path.join(d, 'cube.svg')).read()).replace('%%SECANT%%', open(os.path.join(d, 'secant.svg')).read())
open(os.path.join(d, 'index.html'), 'w').write(t)
