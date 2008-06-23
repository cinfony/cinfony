import sys

if sys.platform[:4] == "java":
    from cdkjython import *
else:
    from cdkjpype import *
