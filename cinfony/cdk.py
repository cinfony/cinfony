#-*. coding: utf-8 -*-
## Copyright (c) 2008-2015, Noel O'Boyle; 2012, Adrià Cereto-Massagué
## All rights reserved.
##
##  This file is part of Cinfony.
##  The contents are covered by the terms of the BSD license
##  which is included in the file LICENSE_BSD.txt.

"""
cdk - A Cinfony module for accessing the CDK from CPython and Jython

Global variables:
  cdk - the underlying CDK Java library (org.openscience.cdk)
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  descs - a list of supported descriptors
  fps - a list of supported fingerprint types
  forcefields - a list of supported forcefields
"""
import sys
import os

if sys.platform[:4] == "java":
    import org.openscience.cdk as cdk
    import java
    import javax

    #Exceptions are handled differently in jpype and jython. We need to wrap them:
    InvalidSmilesException = cdk.exception.InvalidSmilesException
    CDKException = cdk.exception.CDKException
    NullPointerException = java.lang.NullPointerException

else:
    from jpype import *

    if not isJVMStarted():
        _jvm = os.environ['JPYPE_JVM']
        if _jvm[0] == '"': # Remove trailing quotes
            _jvm = _jvm[1:-1]
        _cp = os.environ['CLASSPATH']
        startJVM(_jvm, "-Djava.class.path=" + _cp)

    cdk = JPackage("org").openscience.cdk
    try:
        _testmol = cdk.Molecule()
    except TypeError:
        raise ImportError, "The CDK Jar file cannot be found."

    #Exception wrappers for Jpype
    InvalidSmilesException = JavaException
    CDKException = JavaException
    NullPointerException = JavaException

_chemobjbuilder = cdk.silent.SilentChemObjectBuilder.getInstance()
_aromaticityperceptor = cdk.aromaticity.Aromaticity(cdk.aromaticity.ElectronDonation.daylight(), cdk.graph.Cycles.or(cdk.graph.Cycles.all(), cdk.graph.Cycles.all(6)))

def _getdescdict():
    de = cdk.qsar.DescriptorEngine(cdk.qsar.IMolecularDescriptor, _chemobjbuilder)
    descdict = {}
    for desc in de.getDescriptorInstances():
        spec = desc.getSpecification()
        descclass = de.getDictionaryClass(spec)
        if "proteinDescriptor" not in descclass:
            # Using str() for unicode conversion
            name = str(spec.getSpecificationReference().split("#")[-1])
            descdict[name] = desc
    return descdict

_descdict = _getdescdict()
descs = _descdict.keys()
"""A list of supported descriptors"""
_fingerprinters = {"daylight":cdk.fingerprint.Fingerprinter
                            , "graph":cdk.fingerprint.GraphOnlyFingerprinter
                            , "maccs":cdk.fingerprint.MACCSFingerprinter
                            , "estate":cdk.fingerprint.EStateFingerprinter
                            , "extended":cdk.fingerprint.ExtendedFingerprinter
                            , "hybridization":cdk.fingerprint.HybridizationFingerprinter
                            , "klekota-roth":cdk.fingerprint.KlekotaRothFingerprinter
                            , "pubchem":cdk.fingerprint.PubchemFingerprinter
                            , "substructure":cdk.fingerprint.SubstructureFingerprinter
                            }
fps = _fingerprinters.keys()
"""A list of supported fingerprint types"""
_formats = {'smi': "SMILES" , 'can': "Canonical SMILES", 'sdf': "MDL SDF",
            'mol2': "MOL2", 'mol': "MDL MOL",
            "inchi":"InChI",
            "inchikey":"InChIKey"}
_informats = {'sdf': cdk.io.MDLV2000Reader, 'mol': cdk.io.MDLV2000Reader}
informats = dict([(_x, _formats[_x]) for _x in ['smi', 'sdf', 'mol', 'inchi']])
"""A dictionary of supported input formats"""
_outformats = {'mol': cdk.io.MDLV2000Writer,
               'mol2': cdk.io.Mol2Writer,
               'sdf': cdk.io.SDFWriter}
outformats = dict([(_x, _formats[_x]) for _x in _outformats.keys() + ['can', 'smi', 'inchi', 'inchikey']])
"""A dictionary of supported output formats"""
forcefields = list(cdk.modeling.builder3d.ModelBuilder3D.getInstance(_chemobjbuilder).getFfTypes())
"""A list of supported forcefields"""

_isofact = cdk.config.Isotopes.getInstance()

_bondtypes = {1: cdk.CDKConstants.BONDORDER_SINGLE,
              2: cdk.CDKConstants.BONDORDER_DOUBLE,
              3: cdk.CDKConstants.BONDORDER_TRIPLE}
_revbondtypes = dict([(_y,_x) for (_x,_y) in _bondtypes.iteritems()])

def _intvalue(integer):
    """Paper over some differences between JPype and Jython"""
    # Jython automagically converts Integer to ints
    if type(integer) != type(42): # Is it a Python int?
        integer = integer.intValue()
    return integer

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
    format = format.lower()
    if not os.path.isfile(filename):
        raise IOError, "No such file: '%s'" % filename
    builder = cdk.DefaultChemObjectBuilder.getInstance()
    if format=="sdf":
        return (Molecule(mol) for mol in cdk.io.iterator.IteratingSDFReader(
               java.io.FileInputStream(java.io.File(filename)),
               builder)
               )
    elif format=="smi":
        return (Molecule(mol) for mol in cdk.io.iterator.IteratingSmilesReader(
            java.io.FileInputStream(java.io.File(filename)),
            builder
            ))
    elif format == 'inchi':
        inputfile = open(filename, 'rb')
        return (readstring('inchi', line.rstrip()) for line in inputfile)
    elif format in informats:
        reader = _informats[format](java.io.FileInputStream(java.io.File(filename)))
        chemfile = reader.read(cdk.ChemFile())
        manip = cdk.tools.manipulator.ChemFileManipulator
        return iter(Molecule(manip.getAllAtomContainers(chemfile)[0]),)
    else:
        raise ValueError,"%s is not a recognised CDK format" % format

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
    if format=="smi":
        sp = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
        try:
            ans = sp.parseSmiles(string)
        except InvalidSmilesException, ex:
            if sys.platform[:4] != "java":
                #Jpype exception
                ex = ex.message()
            raise IOError, ex
        return Molecule(ans)
    elif format == 'inchi':
        factory = cdk.inchi.InChIGeneratorFactory.getInstance()
        intostruct = factory.getInChIToStructure(string,cdk.DefaultChemObjectBuilder.getInstance())
        return Molecule(intostruct.getAtomContainer())
    elif format in informats:
        reader = _informats[format](java.io.StringReader(string))
        chemfile = reader.read(cdk.ChemFile())
        manip = cdk.tools.manipulator.ChemFileManipulator
        return Molecule(manip.getAllAtomContainers(chemfile)[0])
    else:
        raise ValueError,"%s is not a recognised CDK format" % format

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
        self.format = format.lower()
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % self.filename
        if not format in outformats:
            raise ValueError,"%s is not a recognised CDK format" % format
        if self.format in ('smi','inchi', 'inchikey'):
            self._outputfile = open(self.filename, "w")
        else:
            self._writer = java.io.FileWriter(java.io.File(self.filename))
            self._molwriter = _outformats[self.format](self._writer)
        self.total = 0 # The total number of molecules written to the file

    def write(self, molecule):
        """Write a molecule to the output file.

        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError, "Outputfile instance is closed."
        if self.format in ('smi','inchi', 'inchikey'):
            self._outputfile.write("%s\n" % molecule.write(format))
        else:
            self._molwriter.write(molecule.Molecule)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.filename = None
        if self.format in ('can', 'smi', 'inchi', 'inchikey'):
            self._outputfile.close()
        else:
            self._molwriter.close()
            self._writer.close()

class Molecule(object):
    """Represent a cdkjpype Molecule.

    Required parameters:
       Molecule -- a CDK Molecule or any type of cinfony Molecule

    Attributes:
       atoms, data, exactmass, formula, molwt, title

    Methods:
       addh(), calcfp(), calcdesc(), draw(), removeh(), write()

    The underlying CDK Molecule can be accessed using the attribute:
       Molecule
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
        _aromaticityperceptor.apply(self.Molecule)

    @property
    def atoms(self): return [Atom(self.Molecule.getAtom(i)) for i in range(self.Molecule.getAtomCount())]
    @property
    def data(self): return MoleculeData(self.Molecule)
    @property
    def formula(self):
        manip = cdk.tools.manipulator.MolecularFormulaManipulator
        mf = manip.getMolecularFormula(self.Molecule)
        return manip.getString(mf) # GetHillString
    @property
    def exactmass(self):
        clone = Molecule(self.Molecule.clone())
        clone.addh()
        manip = cdk.tools.manipulator.MolecularFormulaManipulator
        mf = manip.getMolecularFormula(clone.Molecule)
        return manip.getMajorIsotopeMass(mf)
    @property
    def molwt(self):
        clone = Molecule(self.Molecule.clone())
        clone.addh()
        atommanip = cdk.tools.manipulator.AtomContainerManipulator
        return atommanip.getNaturalExactMass(clone.Molecule)
    def _gettitle(self): return self.Molecule.getProperty(cdk.CDKConstants.TITLE)
    def _settitle(self, val): self.Molecule.setProperty(cdk.CDKConstants.TITLE, val)
    title = property(_gettitle, _settitle)
    @property
    def _exchange(self):
        gt = cdk.geometry.GeometryTools
        if gt.has2DCoordinates(self.Molecule) or gt.has3DCoordinates(self.Molecule):
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
        atommanip = cdk.tools.manipulator.AtomContainerManipulator
        atommanip.convertImplicitToExplicitHydrogens(self.Molecule)

    def removeh(self):
        """Remove hydrogens."""
        atommanip = cdk.tools.manipulator.AtomContainerManipulator
        self.Molecule = atommanip.suppressHydrogens(self.Molecule)

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
        if format not in outformats:
            raise ValueError,"%s is not a recognised CDK format" % format

        if filename is not None and not overwrite and os.path.isfile(filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename

        if format in ("smi", "can"):
            if format == "can":
                sg = cdk.smiles.SmilesGenerator.absolute().aromatic()
            else:
                sg = cdk.smiles.SmilesGenerator.isomeric().aromatic()
            # Set flag or else c1ccccc1 will be written as C1CCCCC1
            smiles = sg.create(self.Molecule)
            if filename:
                output = open(filename, "w")
                print >> output, smiles
                output.close()
                return
            else:
                return smiles
        elif format in ('inchi', 'inchikey'):
            factory = cdk.inchi.InChIGeneratorFactory.getInstance()
            gen = factory.getInChIGenerator(self.Molecule)
            if format == 'inchi':
                return gen.getInchi()
            else:
                return gen.getInchiKey()

        else:
            if filename is None:
                writer = java.io.StringWriter()
            else:
                writer = java.io.FileWriter(java.io.File(filename))
            molwriter = _outformats[format](writer)
            molwriter.write(self.Molecule)
            molwriter.close()
            writer.close()
            if filename == None:
                return str(writer.toString())

    def calcfp(self, fp="daylight"):
        """Calculate a molecular fingerprint.

        Optional parameters:
           fptype -- the fingerprint type (default is "daylight"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        fp = fp.lower()
        if fp in _fingerprinters:
            fingerprinter = _fingerprinters[fp]()
        else:
            raise ValueError, "%s is not a recognised CDK Fingerprint type" % fp
        return Fingerprint(fingerprinter.getBitFingerprint(self.Molecule).asBitSet())

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
            try:
                desc = _descdict[descname]
            except KeyError:
                raise ValueError, "%s is not a recognised CDK descriptor type" % descname
            try:
                value = desc.calculate(self.Molecule).getValue()
                if hasattr(value, "get"): # Instead of array
                    for i in range(value.length()):
                        ans[descname + ".%d" % i] = value.get(i)
                elif hasattr(value, "doubleValue"):
                    ans[descname] = value.doubleValue()
                else:
                    ans[descname] = _intvalue(value)
            except CDKException, ex:
                # Can happen if molecule has no 3D coordinates
                pass
            except NullPointerException, ex:
                # Happens with moment of inertia descriptor
                pass
        return ans

    def draw(self, show=True, filename=None, update=False,
             usecoords=False):
        """Create a 2D depiction of the molecule.

        There is no option to display or write an image file of
        the depiction. For this, you should use the CDK from
        Jython or else the depiction engine of one of the other
        toolkits.

        When using jpype, arguments will be ignored: calling this function is
        equivalent to calling the draw() method of one of the other Cinfony
        modules with parameters:
           show=False, filename=None, update=True, usecoords=False
        """
        if sys.platform[:4] != "java":
            show=False
            filename=None
            update=True
            usecoords=False

        mol = Molecule(self.Molecule.clone())

        if not usecoords:
            # Do the SDG
            sdg = cdk.layout.StructureDiagramGenerator()
            sdg.setUseTemplates(False) # According to JM
            sdg.setMolecule(mol.Molecule)
            sdg.generateCoordinates()
            mol = Molecule(sdg.getMolecule())
            if update:
                for atom, newatom in zip(self.atoms, mol.atoms):
                    coords = newatom.Atom.getPoint2d()
                    atom.Atom.setPoint3d(javax.vecmath.Point3d(
                                         coords.x, coords.y, 0.0))
        else:
           if self.atoms[0].Atom.getPoint2d() is None:
                # Use the 3D coords to set the 2D coords
                for atom, newatom in zip(self.atoms, mol.atoms):
                    coords = atom.Atom.getPoint3d()
                    newatom.Atom.setPoint2d(javax.vecmath.Point2d(
                                    coords.x, coords.y))

        if sys.platform[:4] != "java":
            #We are done in jpype
            return
        mol.removeh()
        canvas = _Canvas(mol.Molecule)

        if filename:
            canvas.writetofile(filename)
        if show:
            canvas.popup()
        else:
            canvas.frame.dispose()

if sys.platform[:4] == "java":
    class _Canvas(javax.swing.JPanel):
        """
        Class used by Molecule.draw() in jython
        """
        def __init__(self, mol):
            self.mol = mol

            self.frame = javax.swing.JFrame()
            generators = []
            generators.append(cdk.renderer.generators.BasicSceneGenerator())
            generators.append(cdk.renderer.generators.standard.StandardGenerator(java.awt.Font("Verdana", java.awt.Font.PLAIN, 18)))
            self.renderer = cdk.renderer.AtomContainerRenderer(generators,
                                            cdk.renderer.font.AWTFontManager())
            model = self.renderer.getRenderer2DModel()
            model.set(cdk.renderer.generators.standard.StandardGenerator.AtomColor, cdk.renderer.color.CDK2DAtomColors())
            drawArea = java.awt.Rectangle(300, 300)
            self.renderer.setup(mol, drawArea)
            image = java.awt.image.BufferedImage(300, 300,
                            java.awt.image.BufferedImage.TYPE_INT_RGB)
            screenSize = java.awt.Dimension(300, 300)
            self.setPreferredSize(screenSize)
            self.setBackground(java.awt.Color.WHITE)
            self.frame.getContentPane().add(self)
            self.frame.pack()
            self.frame.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)

        def paint(self, g):
            javax.swing.JPanel.paint(self, g)
            self.renderer.paint(self.mol, cdk.renderer.visitor.AWTDrawVisitor(g),
                                java.awt.Rectangle(300, 300), True);

        def popup(self):
            self.frame.visible = True

        def writetofile(self, filename):
            img = self.createImage(300, 300)
            g2 = img.getGraphics() # Graphics2D
            g2.setColor(java.awt.Color.WHITE)
            g2.fillRect(0, 0, 300, 300)
            self.paint(g2)
            javax.imageio.ImageIO.write(img, "png", java.io.File(filename))

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
        return cdk.similarity.Tanimoto.calculate(self.fp, other.fp)
    def __getattr__(self, attr):
        if attr == "bits":
            # Create a bits attribute on-the-fly
            bits = []
            idx = self.fp.nextSetBit(0)
            while idx >= 0:
                bits.append(idx)
                idx = self.fp.nextSetBit(idx + 1)
            return bits
        else:
            raise AttributeError, "Fingerprint has no attribute %s" % attr
    def __str__(self):
        return self.fp.toString()

class Atom(object):
    """Represent a cdkjpype Atom.

    Required parameters:
       Atom -- a CDK Atom

    Attributes:
       atomicnum, coords, formalcharge

    The original CDK Atom can be accessed using the attribute:
       Atom
    """

    def __init__(self, Atom):
        self.Atom = Atom

    @property
    def atomicnum(self):
        _isofact.configure(self.Atom)
        return _intvalue(self.Atom.getAtomicNumber())
    @property
    def coords(self):
        coords = self.Atom.point3d
        if not coords:
            coords = self.Atom.point2d
            if not coords:
                return (0., 0., 0.)
        else:
            return (coords.x, coords.y, coords.z)
    @property
    def formalcharge(self):
        _isofact.configure(self.Atom)
        return _intvalue(self.Atom.getFormalCharge())

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
        self.smarts = cdk.smiles.smarts.SMARTSQueryTool(smartspattern, _chemobjbuilder)

    def findall(self, molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.

        Required parameters:
           molecule
        """
        match = self.smarts.matches(molecule.Molecule)
        return list(self.smarts.getUniqueMatchingAtoms())

class MoleculeData(object):
    """Store molecule data in a dictionary-type object

    Required parameters:
      Molecule -- a CDK Molecule

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
        self._mol = Molecule
    def _data(self):
        return self._mol.getProperties()
    def _testforkey(self, key):
        if not key in self:
            raise KeyError, "'%s'" % key
    def keys(self):
        return list(self._data().keySet())
    def values(self):
        return list(self._data().values())
    def items(self):
        return [(k, self[k]) for k in self._data().keySet()]
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self._data())
    def __contains__(self, key):
        return key in self._data()
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

if __name__=="__main__": #pragma: no cover
    mol = readstring("smi", "CC(=O)Cl")
    mol.title = "Noel"
    mol.draw()

    for mol in readfile("sdf", "head.sdf"):
        pass
