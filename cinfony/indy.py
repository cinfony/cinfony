#-*. coding: utf-8 -*-
## Copyright (c) 2011, Noel O'Boyle; 2012, Adrià Cereto-Massagué
## All rights reserved.
##
##  This file is part of Cinfony.
##  The contents are covered by the terms of the GPL v3 license
##  which is included in the file LICENSE_GPLv3.txt.


"""
indy - A Cinfony module for accessing Indigo from CPython, Jython or IronPython

Global variables:
  indigo - the underlying Indigo() object
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  fps - a list of supported fingerprint types
"""

import os
import sys
import tempfile

if sys.platform[:3] == "cli":
    _indigonet = os.environ["INDIGONET"]
    import clr
    clr.AddReference('System.Windows.Forms')
    clr.AddReference('System.Drawing')
    clr.AddReferenceToFileAndPath(_indigonet + "\\indigo-dotnet.dll")
    clr.AddReferenceToFileAndPath(_indigonet + "\\indigo-inchi-dotnet.dll")
    clr.AddReferenceToFileAndPath(_indigonet + "\\indigo-renderer-dotnet.dll")
    from System.Windows.Forms import (
        Application, DockStyle, Form, PictureBox, PictureBoxSizeMode
        )
    from System.Drawing import Image, Size
elif sys.platform[:4] == "java":
    import java, javax

if sys.platform[:3] == "cli" or sys.platform[:4] == "java":
    from com.ggasoftware.indigo import Indigo, IndigoException, IndigoRenderer, IndigoInchi
else:
    from indigo import Indigo, IndigoException
    from indigo_renderer import IndigoRenderer
    from indigo_inchi import IndigoInchi

indigo = Indigo()
indigoInchi = IndigoInchi(indigo)

# PIL and Tkinter
try:
    import Tkinter as tk
    import Image as PIL
    import ImageTk as PILtk
except:
    PILtk = None

fps = ["sim", "sub", "sub-res", "sub-tau", "full"]
"""A list of supported fingerprint types"""

_formats = {'smi': "SMILES", 'can': "Canonical SMILES", "rdf": "MDL RDF file",
            'mol': "MDL MOL file", 'sdf': "MDL SDF file",
            'cml': "Chemical Markup Language",
            'inchi': "InChI", 'inchikey': "InChIKey"}
informats = dict([(_x, _formats[_x]) for _x in ['mol', 'sdf', 'rdf', 'smi',
                                                'cml', 'inchi']])
"""A dictionary of supported input formats"""
outformats = dict([(_x, _formats[_x]) for _x in ['mol', 'sdf', 'smi', 'can',
                                                 'cml', 'inchi', 'inchikey']])
"""A dictionary of supported output formats"""

def readfile(format, filename):
    """Iterate over the molecules in a file.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       filename

    You can access the first molecule in a file using the next() method
    of the iterator:
        mol = readfile("smi", "myfile.smi").next()

    You can make a list of the molecules in a file using:
        mols = list(readfile("smi", "myfile.smi"))

    You can iterate over the molecules in a file as shown in the
    following code snippet:
    >>> atomtotal = 0
    >>> for mol in readfile("sdf", "head.sdf"):
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    if not os.path.isfile(filename):
        raise IOError, "No such file: '%s'" % filename
    format = format.lower()
    # Eagerly evaluate the supplier functions in order to report
    # errors in the format and errors in opening the file.
    # Then switch to an iterator...
    if format=="sdf":
        iterator = indigo.iterateSDFile(filename)
        def sdf_reader():
            for mol in iterator:
                yield Molecule(mol)
        return sdf_reader()
    elif format=="rdf":
        iterator = indigo.iterateRDFile(filename)
        def rdf_reader():
            for mol in iterator:
                yield Molecule(mol)
        return rdf_reader()
    elif format=="mol":
        def mol_reader():
            yield Molecule(indigo.loadMoleculeFromFile(filename))
        return mol_reader()
    elif format=="smi":
        iterator = iterateSmilesFile(filename)
        def smi_reader():
            for mol in iterator:
                yield Molecule(mol)
        return smi_reader()
    elif format=="cml":
        iterator = iterateCMLFile(filename)
        def cml_reader():
            for mol in iterator:
                yield Molecule(mol)
        return cml_reader()

    else:
        raise ValueError, "%s is not a recognised Indigo format" % format

def readstring(format, string):
    """Read in a molecule from a string.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       string

    Example:
    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi", input)
    >>> len(mymol.atoms)
    5
    """
    format = format.lower()
    if format not in informats:
        raise ValueError,"%s is not a recognised Indigo format" % format

    module = indigo if format != "inchi" else indigoInchi
    try:
        mol = module.loadMolecule(string)
    except IndigoException:
        raise IOError, "Failed to convert '%s' to format '%s'" % (
            string, format)

    return Molecule(mol)


class Outputfile(object):
    """Represent a file to which *output* is to be sent.

    Required parameters:
       format - see the outformats variable for a list of available
                output formats
       filename

    Optional parameters:
       overwite -- if the output file already exists, should it
                   be overwritten? (default is False)

    Methods:
       write(molecule)
       close()
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % self.filename
        if self.format in ["sdf", "cml", "rdf", "smi"]:
            self._writer = indigo.writeFile(self.filename)
        else:
            raise ValueError,"%s is not supported for multimolecule output" % format
        self.total = 0 # The total number of molecules written to the file

        if self.format == "cml":
            self._writer.cmlHeader()
        elif self.format == "rdf":
            self._writer.rdfHeader()

    def write(self, molecule):
        """Write a molecule to the output file.

        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError, "Outputfile instance is closed."

        if self.format == "sdf":
            self._writer.sdfAppend(molecule.Mol)
        elif self.format == "rdf":
            self._writer.rdfAppend(molecule.Mol)
        elif self.format == "cml":
            self._writer.cmlAppend(molecule.Mol)
        elif self.format == "smi":
            self._writer.smilesAppend(molecule.Mol)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        if self.format == "cml":
            self._writer.cmlFooter()
        self._writer.close()
        self.filename = None
        del self._writer

class Molecule(object):
    """Represent an Indigo Molecule.

    Required parameter:
       Mol -- an Indigo Mol or any type of cinfony Molecule

    Attributes:
       atoms, data, molwt, title

    Methods:
       addh(), calcfp(), draw(), localopt(), removeh(),
       write()

    The underlying Indigo Molecule can be accessed using the attribute:
       Mol
    """
    _cinfony = True

    def __init__(self, Mol):
        if hasattr(Mol, "_cinfony"):
            a, b = Mol._exchange
            if a == 0:
                molecule = readstring("smi", b)
            else:
                molecule = readstring("mol", b)
            Mol = molecule.Mol

        self.Mol = Mol

    @property
    def atoms(self): return [Atom(atom) for atom in self.Mol.iterateAtoms()]
    @property
    def data(self): return MoleculeData(self.Mol)
    @property
    def formula(self): return self.Mol.grossFormula()
    @property
    def molwt(self): return self.Mol.molecularWeight()
    def _gettitle(self):
        return self.Mol.name()
    def _settitle(self, val): self.Mol.setName(val)
    title = property(_gettitle, _settitle)
    @property
    def _exchange(self):
        if not self.Mol.hasZCoord():
            return (0, self.write("can"))
        else: # If 3D
            return (1, self.write("mol"))

    def addh(self):
        """Add hydrogens."""
        self.Mol.unfoldHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.Mol.foldHydrogens()

    def write(self, format="smi", filename=None, overwrite=False):
        """Write the molecule to a file or return a string.

        Optional parameters:
           format -- see the informats variable for a list of available
                     output formats (default is "smi")
           filename -- default is None
           overwite -- if the output file already exists, should it
                       be overwritten? (default is False)

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.

        To write multiple molecules to the same file you should use
        the Outputfile class.
        """
        format = format.lower()
        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename
        if format=="smi":
            result = self.Mol.smiles()
        elif format=="can":
            result = self.Mol.canonicalSmiles()
        elif format=="mol":
            result = self.Mol.molfile()
        elif format=="inchi":
            result = indigoInchi.getInchi(self.Mol)
        elif format=="inchikey":
            result = indigoInchi.getInchiKey(self.write("inchi"))
        elif format=="cml":
            result = self.Mol.cml()
        elif format=="sdf":
            # No sdf method so use a writeBuffer() as described by Dmitry
            buf = indigo.writeBuffer()
            buf.sdfAppend(self.Mol)
            result = buf.toString()
        else:
            raise ValueError,"%s is not a recognised Indigo format" % format
        if filename:
            output = open(filename, "w")
            output.write(result)
            output.close()
        else:
            return result

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.

        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        return iter(self.atoms)

    def __str__(self):
        return self.write()

##    def calcdesc(self, descnames=[]):
##        """Calculate descriptor values.
##
##        Optional parameter:
##           descnames -- a list of names of descriptors
##
##        If descnames is not specified, all available descriptors are
##        calculated. See the descs variable for a list of available
##        descriptors.
##        """
##        if not descnames:
##            descnames = descs
##        ans = {}
##        for descname in descnames:
##            try:
##                desc = descDict[descname]
##            except KeyError:
##                raise ValueError, "%s is not a recognised RDKit descriptor type" % descname
##            ans[descname] = desc(self.Mol)
##        return ans

    def calcfp(self, fptype="sim"):
        """Calculate a molecular fingerprint.

        Optional parameters:
           fptype -- the fingerprint type (default is "sim"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        fptype = fptype.lower()
        if fptype in ["sim", "sub", "sub-res", "sub-tau", "full"]:
            fp = Fingerprint(self.Mol.fingerprint(fptype))
        else:
            raise ValueError, "%s is not a recognised Indigo Fingerprint type" % fptype
        return fp

    def draw(self, show=True, filename=None, update=False, usecoords=False):
        """Create a 2D depiction of the molecule.

        Optional parameters:
          show -- display on screen (default is True)
          filename -- write to file (default is None)
          update -- update the coordinates of the atoms to those
                    determined by the structure diagram generator
                    (default is False)
          usecoords -- don't calculate 2D coordinates, just use
                       the current coordinates (default is False)

        Tkinter and Python Imaging Library are required for image display.
        """
        if update:
            mol = self.Mol
        else:
            mol = self.Mol.clone()
        if not usecoords:
            mol.layout()
        if show or filename:
            renderer = IndigoRenderer(indigo)
            indigo.setOption("render-output-format", "png")
            indigo.setOption("render-margins", 10, 10)
            indigo.setOption("render-coloring", "True")
            indigo.setOption("render-image-size", 300, 300)
            indigo.setOption("render-background-color", "1.0, 1.0, 1.0")
            if self.title:
                indigo.setOption("render-comment", self.title)
            if filename:
                filedes = None
            else:
                filedes, filename = tempfile.mkstemp()

            renderer.renderToFile(mol, filename)

            if show:
                if sys.platform[:4] == "java":
                    image = javax.imageio.ImageIO.read(java.io.File(filename))
                    frame = javax.swing.JFrame(visible=1)
                    frame.getContentPane().add(javax.swing.JLabel(javax.swing.ImageIcon(image)))
                    frame.setSize(300,300)
                    frame.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
                    frame.show()

                elif sys.platform[:3] == "cli":
                    if filedes:
                        errormessage = ("It is only possible to show the molecule if you "
                                        "provide a filename. The reason for this is that I kept "
                                        "having problems when using temporary files.")
                        raise RuntimeError(errormessage)
                    form = Form()
                    form.ClientSize = Size(300, 300)
                    form.Text = self.title
                    image = Image.FromFile(filename)
                    box = PictureBox()
                    box.SizeMode = PictureBoxSizeMode.StretchImage
                    box.Image = image
                    box.Dock = DockStyle.Fill
                    form.Controls.Add(box)
                    form.Show()
                    Application.Run(form)

                else:
                    if not PILtk:
                        errormessage = ("Tkinter or Python Imaging "
                                        "Library not found, but is required for image "
                                        "display. See installation instructions for "
                                        "more information.")
                        raise ImportError, errormessage

                    root = tk.Tk()
                    root.title((hasattr(self, "title") and self.title)
                               or self.__str__().rstrip())
                    frame = tk.Frame(root, colormap="new", visual='truecolor').pack()
                    image = PIL.open(filename)
                    imagedata = PILtk.PhotoImage(image)
                    label = tk.Label(frame, image=imagedata).pack()
                    quitbutton = tk.Button(root, text="Close", command=root.destroy).pack(fill=tk.X)
                    root.mainloop()


            if filedes:
                os.close(filedes)
                os.remove(filename)

class Atom(object):
    """Represent an Indigo Atom.

    Required parameters:
       Atom -- an Indigo Atom

    Attributes:
        atomicnum, coords, formalcharge

    The original Indigo Atom can be accessed using the attribute:
       Atom
    """

    def __init__(self, Atom):
        self.Atom = Atom
    @property
    def atomicnum(self): return self.Atom.atomicNumber()
    @property
    def coords(self):
        return tuple(self.Atom.xyz())
    @property
    def formalcharge(self): return self.Atom.charge()

    def __str__(self):
        if hasattr(self, "coords"):
            return "Atom: %d (%.2f %.2f %.2f)" % (self.atomicnum, self.coords[0],
                                                    self.coords[1], self.coords[2])
        else:
            return "Atom: %d (no coords)" % (self.atomicnum)

class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern

    Methods:
       findall(molecule)

    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol)
    [(0, 1), (3, 4), (5, 6)]

    The numbers returned are the indices (starting from 0) of the atoms
    that match the SMARTS pattern. In this case, there are three matches
    for each of the three ethyl groups in the molecule.
    """
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        try:
            self.smarts = indigo.loadSmarts(smartspattern)
        except IndigoException:
            raise IOError, "Invalid SMARTS pattern."

    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.

        Required parameters:
           molecule
        """
        matcher = indigo.substructureMatcher(molecule.Mol)
        matches = list(matcher.iterateMatches(self.smarts))
        ans = []
        for match in matches:
            a = []
            for queryatom in self.smarts.iterateAtoms():
                a.append(match.mapAtom(queryatom).index())
            ans.append(tuple(a))
        return ans

class MoleculeData(object):
    """Store molecule data in a dictionary-type object

    Required parameters:
      Mol -- an Indigo Mol

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying Mol.

    Example:
    >>> mol = readfile("sdf", 'head.sdf').next()
    >>> data = mol.data
    >>> print data
    {'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}
    >>> print len(data), data.keys(), data.has_key("NSC")
    2 ['Comment', 'NSC'] True
    >>> print data['Comment']
    CORINA 2.61 0041  25.10.2001
    >>> data['Comment'] = 'This is a new comment'
    >>> for k,v in data.iteritems():
    ...    print k, "-->", v
    Comment --> This is a new comment
    NSC --> 1
    >>> del data['NSC']
    >>> print len(data), data.keys(), data.has_key("NSC")
    1 ['Comment'] False
    """
    def __init__(self, Mol):
        self._mol = Mol
    def _testforkey(self, key):
        if not self._mol.hasProperty(key):
            raise KeyError, "'%s'" % key
    def keys(self):
        return [prop.name() for prop in self._mol.iterateProperties()]
    def values(self):
        return [prop.rawData() for prop in self._mol.iterateProperties()]
    def items(self):
        return [(prop.name(), prop.rawData())
                for prop in self._mol.iterateProperties()]
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self.keys())
    def __contains__(self, key):
        return self._mol.hasProperty(key)
    def __delitem__(self, key):
        self._testforkey(key)
        self._mol.removeProperty(key)
    def clear(self):
        for key in self:
            del self[key]
    def has_key(self, key):
        return key in self
    def update(self, dictionary):
        for k, v in dictionary.iteritems():
            self[k] = v
    def __getitem__(self, key):
        self._testforkey(key)
        return self._mol.getProperty(key)
    def __setitem__(self, key, value):
        self._mol.setProperty(key, str(value))
    def __repr__(self):
        return dict(self.iteritems()).__repr__()

class Fingerprint(object):
    """A Molecular Fingerprint.

    Required parameters:
       fingerprint -- a vector calculated by one of the fingerprint methods

    Attributes:
       fp -- the underlying fingerprint object
       bits -- a list of bits set in the Fingerprint

    Methods:
       The "|" operator can be used to calculate the Tanimoto coeff. For example,
       given two Fingerprints 'a', and 'b', the Tanimoto coefficient is given by:
          tanimoto = a | b
    """
    def __init__(self, fingerprint):
        self.fp = fingerprint
    def __or__(self, other):
        return indigo.similarity(self.fp, other.fp, "tanimoto")
    def _buffer_to_int(self):
        stringrep = self.fp.toString()
        return [int(stringrep[i:i+1]) for i in range(0, len(stringrep), 1)]
    @property
    def bits(self):
        return _findbits(self._buffer_to_int(), 8)
    def __str__(self):
        return str(self._buffer_to_int())

def _toint(string):
    """
    Some bits sometimes are a character. I haven't found what do they mean,
    but they break cinfony fingerprints unless taken care of. This functions is just for that.
    """
    if string.isdigit():
        return int(string)
    else:
        return 0

def _findbits(fp, bitsperint):
    """Find which bits are set in a list/vector.

    This function is used by the Fingerprint class.

    >>> _findbits([13, 71], 8)
    [1, 3, 4, 9, 10, 11, 15]
    """
    ans = []
    start = 1
    for x in fp:
        i = start
        while x > 0:
            if x % 2:
                ans.append(i)
            x >>= 1
            i += 1
        start += bitsperint
    return ans


def _compressbits(bitvector, wordsize=32):
    """Compress binary vector into vector of long ints.

    This function is used by the Fingerprint class.

    >>> _compressbits([0, 1, 0, 0, 0, 1], 2)
    [2, 0, 2]
    """
    ans = []
    for start in range(0, len(bitvector), wordsize):
        compressed = 0
        for i in range(wordsize):
            if i + start < len(bitvector) and bitvector[i + start]:
                compressed += 2**i
        ans.append(compressed)

    return ans

if __name__=="__main__": #pragma: no cover
    import doctest
    doctest.testmod()

