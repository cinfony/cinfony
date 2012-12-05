#-*. coding: utf-8 -*-
## Copyright (c) 2012, Adrià Cereto-Massagué, Noel O'Boyle
## All rights reserved.
##
##  This file is part of Cinfony.
##  The contents are covered by the terms of the BSD license
##  which is included in the file LICENSE_BSD.txt.

"""
jchem - A Cinfony module for accessing ChemAxon's JChem from CPython and Jython

Global variables:
  chemaxon - the underlying JChem Java library
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  descs - a list of supported descriptors
  fps - a list of supported fingerprint types
  forcefields - a list of supported forcefields
"""
import sys
import os
from glob import glob

if sys.platform[:4] == "java":
    classpath = []
    if 'JCHEMDIR' in os.environ:
        assert os.path.isdir(os.path.join(os.environ['JCHEMDIR'], 'lib'))
        for jar in glob(os.path.join(os.path.join(os.environ['JCHEMDIR'],'lib'), '*.jar')):
            classpath.append(jar)

if sys.platform[:4] == "java" or sys.platform[:3] == "cli":
    import sys
    sys.path = classpath + sys.path
    import java, javax
    import chemaxon
    from chemaxon.util import MolHandler
    #Exceptions are handled differently in jpype and jython. We need to wrap them:
    MolExportException = chemaxon.marvin.io.MolExportException
    MolFormatException = chemaxon.formats.MolFormatException
else:
    from jpype import *

    if not isJVMStarted():
        _jvm = os.environ['JPYPE_JVM']
        if _jvm[0] == '"': # Remove trailing quotes
            _jvm = _jvm[1:-1]
        _cp = os.pathsep.join(os.environ.get('CLASSPATH', '').split(os.pathsep))
        startJVM(_jvm, "-Djava.class.path=" + _cp)

    chemaxon = JPackage("chemaxon")
    MolHandler = chemaxon.util.MolHandler
    try:
        _testmol = MolHandler()
    except TypeError:
        raise ImportError, "jchem.jar file cannot be found."

    # Exception wrappers for JPype
    MolExportException = JavaException
    MolFormatException = JavaException

_descset = set(['HAcc', 'HDon', 'Heavy', 'LogD', 'LogP', 'Mass', 'TPSA'])
_descset.update(dir(chemaxon.descriptors.scalars))
descs = [cls for cls in _descset if hasattr(getattr(chemaxon.descriptors.scalars, cls),'generate') and cls != 'LogD'] + ['RotatableBondsCount']
"""A list of supported descriptors"""
fps = ['ecfp']
"""A list of supported fingerprint types"""
forcefields = ["mmff94"]
"""A list of supported forcefields"""

informats = {
    'smi': "SMILES"
    ,'cxsmi': "ChemAxon exntended SMILES"
    ,'mol': "MDL MOL"
    ,'sdf': "MDL SDF"
    ,'inchi': "InChI"
    ,'cml': "Chemical Markup Language"
    , 'mrv':'Marvin Documents'
    , 'skc':'ISIS/Draw sketch file'
    , 'cdx':'ChemDraw sketch file'
    , 'cdxml':'ChemDraw sketch file'
    , "name":"Common name"
    , "peptide":"Aminoacid sequence"
    , "sybyl":"Tripos SYBYL"
    , "pdb":"PDB"
    , "xyz":"XYZ"
    , 'cube':'Gaussian cube'
    , 'gout':'Gaussian output format'
    }
"""A dictionary of supported input formats"""

outformats = {
    'smi': "SMILES"
    ,'cxsmi': "ChemAxon exntended SMILES"
    ,'mol': "MDL MOL"
    ,'sdf': "MDL SDF"
    ,'inchi': "InChI"
    ,'inchikey': "InChIKey"
    ,'cml': "CML"
    , 'mrv':'Marvin Documents'
    , 'skc':'ISIS/Draw sketch file'
    , 'cdx':'ChemDraw sketch file'
    , 'cdxml':'ChemDraw sketch file'
    , "name":"Common name"
    , "peptide":"Aminoacid sequence"
    , "sybyl":"Tripos SYBYL"
    , "pdb":"PDB"
    , "xyz":"XYZ"
    , 'cube':'Gaussian cube'
    , 'gjf':'Gaussian input format'
    }
"""A dictionary of supported output formats"""

def readfile(format, filename):
    """Iterate over the molecules in a file.

    Required parameters:
       format - Ignored, but needed for compatibility with other cinfony
                modules and also good for readability
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
    if not format in outformats:
        raise ValueError("%s is not a recognised JChem format" % format)
    try:
        mi = chemaxon.formats.MolImporter(filename)
        mol = mi.read()
        while mol:
            mol.aromatize()
            yield Molecule(mol)
            mol = mi.read()
    except chemaxon.formats.MolFormatException:
        raise ValueError("%s is not a recognised JChem format" % format)

def readstring(format, string):
    """Read in a molecule from a string.

    Required parameters:
       format - Ignored, but needed for compatibility with other cinfony
                modules and also good for readability
       string

    Example:
    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi", input)
    >>> len(mymol.atoms)
    5
    """
    format = format.lower()
    if format not in informats:
        raise ValueError("%s is not a recognised JChem format" % format)
    try:
        mh = MolHandler(string)
        return Molecule(mh.molecule)
    except MolFormatException, ex:
        if sys.platform[:4] != "java":
            #Jpype exception
            ex = ex.message()
            raise IOError, ex
        else:
            raise IOError("Problem reading the supplied string")

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
        if ':' in format:
            format,  options = format.split(':')
            if options:
                options = ':' + options
        else:
            options = ''
        self.format = format.lower()
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % self.filename
        if format in ("smi", 'cxsmi'):
            if not options:
                options = ':a-H'
            out = chemaxon.formats.MolExporter.exportToFormat(self.Molecule,format +'les:a-H')
        try:
            self._writer = chemaxon.formats.MolExporter(filename, format + options)
        except MolExportException,  e:
            raise ValueError(e)
        self.total = 0 # The total number of molecules written to the file

    def write(self, molecule):
        """Write a molecule to the output file.

        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError, "Outputfile instance is closed."
        self._writer.write(molecule.Molecule)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.filename = None
        self._writer.close()

class Molecule(object):
    """Represent a JChem Molecule.

    Required parameters:
       Molecule -- a JChem Molecule or any type of cinfony Molecule

    Attributes:
       atoms, data, exactmass, formula, molwt, title

    Methods:
       addh(), calcfp(), calcdesc(), draw(), removeh(), write()

    The underlying JChem Molecule can be accessed using the attribute:
       Molecule
    The associated JChem MolHandler can be accessed using the attribute:
       MolHandler
    """
    _cinfony = True

    def __init__(self, Molecule):

        if hasattr(Molecule, "_cinfony"):
            a, b = Molecule._exchange
            if a == 0:
                mol = readstring("smi", b)
            else:
                mol = readstring("sdf", b)
            Molecule = mol.Molecule

        self.Molecule = Molecule
        self.MolHandler = chemaxon.util.MolHandler(self.Molecule)
        self.MolHandler.aromatize()

    @property
    def atoms(self): return [Atom(atom) for atom in self.Molecule.atomArray]
    @property
    def data(self): return MoleculeData(self)
    @property
    def formula(self): return self.MolHandler.calcMolFormula()
    @property
    def exactmass(self):
        return self.MolHandler.calcMolWeightInDouble()
    @property
    def molwt(self):
        return self.MolHandler.calcMolWeight()
    def _gettitle(self): return self.Molecule.getName()
    def _settitle(self, val): self.Molecule.setName(val)
    title = property(_gettitle, _settitle)
    @property
    def _exchange(self):
        if self.Molecule.dim > 1:
            return (1, self.write("mol"))
        else:
            return (0, self.write("smi"))

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.

        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        return iter(self.atoms)

    def __str__(self):
        return self.write()

    def addh(self):
        """Add hydrogens."""
        self.MolHandler.addHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.MolHandler.removeHydrogens()

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
        if ':' in format:
            format,  options = format.split(':')
            if options:
                options = ':' + options
        else:
            options = ''
        format = format.lower()
        if format not in outformats:
            raise ValueError("%s is not a recognised format" % format)

        if filename is not None and not overwrite and os.path.isfile(filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename

        if format in ("smi", 'cxsmi'):
            if not options:
                options = ':a-H'
            out = chemaxon.formats.MolExporter.exportToFormat(self.Molecule,format +'les' + options)
        elif format == 'inchikey':
            out = chemaxon.formats.MolExporter.exportToFormat(self.Molecule,'inchikey').replace('InChIKey=', '')
        else:
            out = chemaxon.formats.MolExporter.exportToFormat(self.Molecule,format + options)
            if format == 'inchi':
                out = out.split('AuxInfo=')[0]
        if filename:
            output = open(filename, "w")
            print >> output, out
            output.close()
            return
        else:
            return out


    def calcfp(self, fp="ecfp"):
        """Calculate a molecular fingerprint.

        Optional parameters:
           fptype -- the fingerprint type (default is "daylight"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        fp = fp.lower()
        if fp in fps:
            if fp == 'ecfp':
                fp = chemaxon.descriptors.ECFP(ECFPConfiguration)
                fp.generate(self.Molecule)
        else:
            raise ValueError, "%s is not a recognised fingerprint type" % fp
        return Fingerprint(fp)

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, all available descriptors are
        calculated. See the descs variable for a list of available
        descriptors.
        """
        if not descnames:
            descnames = descs
        ans = {}
        for descname in descnames:
            if descname not in descs:
                raise ValueError, "%s is not a recognised descriptor type" % descname
            if descname == 'RotatableBondsCount':
                ta = chemaxon.calculations.TopologyAnalyser()
                ta.setMolecule(self.Molecule)
                ans[descname] = ta.rotatableBondCount()
            else:
                desc = getattr(chemaxon.descriptors.scalars, descname)('')
                desc.generate(self.Molecule)
                ans[descname] = desc.toFloatArray()[0]
        return ans

    def make3D(self):
        """Generate 3D coordinates.

        Hydrogens are added, and a low energy conformer is found
        using the MMFF94 forcefield.
        """
        self.addh()
        cp = chemaxon.marvin.calculations.ConformerPlugin()
        cp.setMolecule(self.Molecule)
        cp.setLowestEnergyConformerCalculation(True)
        cp.setMMFF94Optimization(True)
        success = cp.run()
        optmol = cp.getMMFF94OptimizedStrucutre()
        self.Molecule = optmol
        self.MolHandler = chemaxon.util.MolHandler(self.Molecule)
        self.MolHandler.aromatize()

    def draw(self, show=True, filename=None, update=False,
             usecoords=False):
        """Create a 2D depiction of the molecule.
        """
        if not usecoords:
            molecule = self.Molecule.clone()
            molecule.setDim(0)
        else:
            molecule = self.Molecule
        if update:
            myMolecule = readstring("mol", Molecule(molecule).write("mol"))
            self.Molecule = myMolecule.Molecule
            self.MolHandler = myMolecule.MolHandler
        bytearray = chemaxon.formats.MolExporter.exportToBinFormat(molecule, 'png')
        if filename:
            of = java.io.FileOutputStream(filename)
            of.write(bytearray)
            of.close()
        if show:
            source = java.io.ByteArrayInputStream(bytearray)
            reader = javax.imageio.ImageIO.getImageReadersByFormatName('png').next()
            iis = javax.imageio.ImageIO.createImageInputStream(source)
            reader.setInput(iis, True)
            param = reader.getDefaultReadParam()
            image = reader.read(0, param)
            frame = javax.swing.JFrame()
            imageIcon = javax.swing.ImageIcon(image)
            label = javax.swing.JLabel()
            label.setIcon(imageIcon)
            frame.getContentPane().add(label, java.awt.BorderLayout.CENTER)
            frame.pack()
            frame.setVisible(True)
            frame.show()


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
        return 1 - self.fp.getTanimoto(other.fp)
    def __getattr__(self, attr):
        if attr == "bits":
            # Create a bits attribute on-the-fly
            bs = self.fp.toBitSet()
            bits = [-1]
            while True:
                setbit = bs.nextSetBit(bits[-1] + 1)
                if setbit == -1:
                    break
                bits.append(setbit)
            return bits[1:] # Leave out the initial '-1'
        else:
            raise AttributeError, "Fingerprint has no attribute %s" % attr
    def __str__(self):
        return ", ".join([str(x) for x in self.fp.toIntArray()])

class Atom(object):
    """Represent an Atom.

    Required parameters:
       Atom -- a JChem Atom

    Attributes:
       atomicnum, coords, formalcharge

    The original JChem Atom can be accessed using the attribute:
       Atom
    """

    def __init__(self, Atom):
        self.Atom = Atom

    @property
    def atomicnum(self): return self.Atom.getAtno()
    @property
    def coords(self):
            return (self.Atom.x, self.Atom.y, self.Atom.z)
    @property
    def formalcharge(self):
        return self.Atom.charge

    def __str__(self):
        c = self.coords
        return "Atom: %d (%.2f %.2f %.2f)" % (self.atomicnum, c[0], c[1], c[2])

class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern

    Methods:
       findall()

    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol)
    [(1, 2), (4, 5), (6, 7)]
    """
    def __init__(self, smartspattern):
        """Initialise with a SMARTS pattern."""
        self.search = chemaxon.sss.search.MolSearch()
        smarts = MolHandler(smartspattern)
        smarts.setQueryMode(True)
        smarts.aromatize()
        self.search.setQuery(smarts.molecule)

    def findall(self, molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.

        Required parameters:
           molecule
        """
        self.search.setTarget(molecule.Molecule)
        match = self.search.findAll()
        result = []
        for i in xrange(len(match)):
            result.append(tuple([n+1 for n in match[i]]))
        return result

class MoleculeData(object):
    """Store molecule data in a dictionary-type object

    Required parameters:
      Molecule -- a JChem Molecule

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying Molecule.

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
    def __init__(self, Molecule):
        self._data = Molecule.Molecule.properties()
    def _testforkey(self, key):
        if not key in self:
            raise KeyError, "'%s'" % key
    def keys(self):
        return list(self._data.keys)
    def values(self):
        return [self[k] for k in self._data.keys]
    def items(self):
        return [(k, self[k]) for k in self._data.keys]
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self._data.keys)
    def __contains__(self, key):
        return key in self.keys()
    def __delitem__(self, key):
        self._testforkey(key)
        self._data.setString(key, None)
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
        return self._data.get(key).propValue
    def __setitem__(self, key, value):
        self._data.setString(key, str(value))
    def __repr__(self):
        return dict(self.iteritems()).__repr__()

ECFPConfiguration = """<?xml version="1.0" encoding="UTF-8"?>
<ECFPConfiguration Version="0.1">

    <Parameters Length="1024" Diameter="4" Counts="no"/>

    <IdentifierConfiguration>
        <!-- Default atom properties (switched on by Value=1) -->
        <Property Name="AtomicNumber" Value="1"/>
        <Property Name="HeavyNeighborCount" Value="1"/>
        <Property Name="HCount" Value="1"/>
        <Property Name="FormalCharge" Value="1"/>
        <Property Name="IsRingAtom" Value="1"/>

        <!-- Other built-in atom properties (switched off by Value=0) -->
        <Property Name="ConnectionCount" Value="0"/>
        <Property Name="Valence" Value="0"/>
        <Property Name="Mass" Value="0"/>
        <Property Name="MassNumber" Value="0"/>
        <Property Name="HasAromaticBond" Value="0"/>
        <Property Name="IsTerminalAtom" Value="0"/>
        <Property Name="IsStereoAtom" Value="0"/>
    </IdentifierConfiguration>

    <StandardizerConfiguration Version="0.1">
        <Actions>
            <Action ID="aromatize" Act="aromatize"/>
            <RemoveExplicitH ID="RemoveExplicitH" Groups="target"/>
        </Actions>
    </StandardizerConfiguration>

    <ScreeningConfiguration>
        <ParametrizedMetrics>
            <ParametrizedMetric Name="Tanimoto" ActiveFamily="Generic" Metric="Tanimoto" Threshold="0.5"/>
            <ParametrizedMetric Name="Euclidean" ActiveFamily="Generic" Metric="Euclidean" Threshold="10"/>
        </ParametrizedMetrics>
    </ScreeningConfiguration>

</ECFPConfiguration>
"""

if __name__=="__main__": #pragma: no cover
    mol = readstring("smi", "CC(=O)Cl")
    mol.title = u"Adrià"
    mol.draw()

    for mol in readfile("sdf", "head.sdf"):
        pass
