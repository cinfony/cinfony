"""
obabel - A Cinfony module for accessing OpenBabel

obabel can be used from all of CPython, Jython and IronPython. It imports the
appropriate Cinfony module (either pybel, jybel or ironable) depending on the
Python implementation.
"""
import sys

if sys.platform[:4] == "java":
    from jybel import *
elif sys.platform[:3] == "cli":
    from ironable import *
else:
    from pybel import *
