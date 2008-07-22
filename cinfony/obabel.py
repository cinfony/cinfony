import sys

if sys.platform[:4] == "java":
    from jybel import *
else:
    from pybel import *
