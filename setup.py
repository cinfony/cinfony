"""cinfony: a common API to several cheminformatics toolkits

cinfony is a Python library that provides a common API to several
open source cheminformatics toolkits.
"""

doclines = __doc__.split("\n")

# My code is BSD
# Open Babel is GPL v2
# OPSIN is Artistic License v2.0 (not viral)
# CDK is LGPL (not viral)
# RDKit is BSD
# Indigo is GPL v3

# Chosen from http://www.python.org/pypi?:action=list_classifiers
classifiers = """\
Development Status :: 5 - Production/Stable
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
License :: OSI Approved :: GNU General Public License (GPL)
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Chemistry
Topic :: Software Development :: Libraries :: Python Modules
"""

from distutils.core import setup
setup(
    name = "cinfony",
    version = "1.1",
    url = "http://cinfony.googlecode.com",
    author = "Noel O'Boyle",
    author_email = "baoilleach@gmail.com",
    maintainer = "Noel O'Boyle",
    maintainer_email = "baoilleach@gmail.com",
    license = "BSD",
    description = doclines[0],
    long_description = "\n".join(doclines[2:]),      
    classifiers = filter(None, classifiers.split("\n")),
    platforms = ["Any."],
    packages = ['cinfony'],
    )

