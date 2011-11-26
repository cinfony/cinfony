## Copyright (c) 2008-2011, Noel O'Boyle
## All rights reserved.
##
##  This file is part of Cinfony.
##  The contents are covered by the terms of the BSD license
##  which is included in the file LICENSE_BSD.txt.

"""
cdk - A Cinfony module for accessing the CDK

cdk can be used from both CPython and Jython. It imports the appropriate
Cinfony module, either cdkjpype and cdkjython, depending on the Python
implementation.
"""
import sys

if sys.platform[:4] == "java":
    from cdkjython import *
else:
    from cdkjpype import *
