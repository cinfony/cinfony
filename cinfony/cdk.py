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
