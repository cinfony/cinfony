import os

try:
    test = os.write
    from cdkjpype import *
except AttributeError:
    from cdkjython import *
