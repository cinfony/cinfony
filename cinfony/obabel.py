"""
obabel - A Cinfony module for accessing OpenBabel

obabel can be used from both CPython and Jython. It imports the appropriate
Cinfony module, either jybel or pybel, depending on the Python implementation.
"""
import sys

if sys.platform[:4] == "java":
    from jybel import *
else:
    from pybel import *
