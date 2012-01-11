## Copyright (c) 2008-2011, Noel O'Boyle
## All rights reserved.
##
##  This file is part of Cinfony.
##  The contents are covered by the terms of the GPL v2 license
##  which is included in the file LICENSE_GPLv2.txt.

"""
obabel - A Cinfony module for accessing OpenBabel

obabel can be used from all of CPython, Jython and IronPython. It imports the
appropriate Cinfony module (either pybel, jybel or ironable) depending on the
Python implementation.
"""
import sys

if sys.platform[:3] == "cli":
    from ironable import *
else:
    from pybel import *
