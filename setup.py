"""cinfony: a common API to several cheminformatics toolkits

cinfony is a Python library that provides a common API to several
open source cheminformatics toolkits.
"""

doclines = __doc__.split("\n")

# Chosen from http://www.python.org/pypi?:action=list_classifiers
classifiers = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Chemistry
Topic :: Software Development :: Libraries :: Python Modules
"""

from distutils.core import setup
setup(
    name = "cinfony",
    version = "0.8.1",
    url = "http://code.googlecode.com/p/cinfony",
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

