from __future__ import generators

import os
import urllib
import tempfile
import StringIO

import org.openscience.cdk as cdk
import java

def _getdescdict():
    de = cdk.qsar.DescriptorEngine(cdk.qsar.DescriptorEngine.MOLECULAR)
    descdict = {}
    for desc in de.getDescriptorInstances():
        spec = desc.getSpecification()
        descclass = de.getDictionaryClass(spec)
        if "proteinDescriptor" not in descclass:
            name = spec.getSpecificationReference().split("#")[-1]
            descdict[name] = desc
    return descdict
_descdict = descdict = _getdescdict()
descriptors = _descdict.keys()

_informats = {}
informats = ['smi' ,'sdf']
_outformats = {'mol': cdk.io.MDLWriter,
               'mol2': cdk.io.Mol2Writer,
               'smi': cdk.io.SMILESWriter,
               'sdf': cdk.io.MDLWriter} # FIXME: Handle the one molecule case
outformats = ['smi'] + _outformats.keys()
forcefields = list(cdk.modeling.builder3d.ModelBuilder3D.getInstance().getFfTypes())

#from org.openscience.cdk.io.iterator import IteratingMDLReader
#from org.openscience.cdk.io import ReaderFactory, WriterFactory, SMILESWriter
#from org.openscience.cdk.smiles import SmilesParser, SmilesGenerator
#from org.openscience.cdk.smiles.smarts import SMARTSParser
#from org.openscience.cdk.isomorphism import UniversalIsomorphismTester

##_informats = dict([(y.formatName, y) for y in ReaderFactory().formats])
##_outformats = dict([(y.formatName, y) for y in WriterFactory().findChemFormats(0)])
##_formatlookup = {
##    "ab": "ABINIT",
##    "ace": "Aces2",
##    "adf": "ADF",
##    "alc": "Alchemy",
##    "bgf": "MSI BGF",
##    "bs": "Ball and Stick",
##    "cac":   "CAChe MolStruct",
##    "cache": "CAChe MolStruct",
##    "cdk": "CDK Source Code",
##    "ctx": "CTX", # Poor description
##    "caccrt": "Cacao Cartesian",
##    "cacint": "Cacao Internal",
##    "c3d1": "Chem3D Cartesian 1",
##    "c3d2": "Chem3D Cartesian 2",
##    "cdx": "ChemDraw eXchange file",
##    "cml": "Chemical Markup Language",
##    "cmlrss": "CML enriched RSS",
##    "crk2d": "Chemical Resource Kit 2D",
##    "crk3d": "Chemical Resource Kit 3D",
##    "cht": "Chemtool",
##    "cry": "CrystClust",
##    "cif": "Crystallographic Interchange Format",
##    "dal": "Dalton",
##    "dmol": "DMol3",
##    "dock": ock 5 Box",
##    "fh": "Fenske-Hall Z-Matrix",
##    "fpt": "Fingerprint",
##    "gam":    "GAMESS log file",
##    "gamout": "GAMESS log file",
##    "gr96": "GROMOS96",
##    "gin": "Gaussian Input",
##    "g90": "Gaussian90",
##    "g92": "Gaussian92",
##    "g94": "Gaussian94",
##    "g95": "Gaussian95",
##    "g98": "Gaussian98",
##    "g03": "Gaussian 2003",
##    "gpr": "Ghemical Quantum/Molecular Mechanics Model",
##    "gsp": 'Ghemical Simplified Protein Model',
##    "hin": "HyperChem HIN",
##    "inchixml": "IUPAC-NIST Chemical Identifier (XML)",
##    "inchi": "IUPAC-NIST Chemical Identifier (Plain Text)",
##    "jag": "Jaguar",
##    "jme": "JME", # Poor description
##    "mol": "MDL Molfile",
##    "mol2000": "MDL Mol/SDF V2000",
##    "mol3000": "MDL Mol/SDF V3000",
##    "rxn": "MDL Reaction format",
##    "rxn3000": "MDL RXN V3000",
##    "sdf": "MDL Structure-data file",
##    "sd":  "MDL Structure-data file",
##    "macie": "MACiE",
##    "mmd": "MacroModel",
##    "mmod": "MacroModel",
##    "mpqc": "Massively Parallel Quantum Chemistry Program",
##    "mol2": "Mol2 (Sybyl)",
##    "mop97": "MOPAC 97",
##    "mop93": "MOPAC 93",
##    "mop7": "MOPAC7",
##    "mop02": "MOPAC 2002",
##    "nw": "NWChem",
##    "pcm": "PCModel",
##    "pqs": "Parallel Quantum Solutions",
##    "pmp": "PolyMorph Predictor (Cerius)",
##    "pdb": "Protein Brookhave Database (PDB)",
##    "pdbml": "Protein Data Bank Markup Language (PDBML)",
##    "ent": "Protein Brookhave Database (PDB)", # Typo
##    "pc": "PubChem",
##    "pcasn": "PubChem Compound ASN",
##    "qc": "Q-Chem",
##    "raw": "Raw Copy",
##    "smi": "SMILES",
##    "fix": "SMILES FIX",
##    "res": "ShelXL",
##    "ins": "ShelXL",
##    "sma": "SMARTS",
##    "spartan": "Spartan Quantum Mechanics Program",
##    "mpd": "Sybyl descriptor",
##    "txyz": "Tinker XYZ",
##    "tmm2": "Tinker MM2", # OpenBabel calls this txyz
##    "tmol": "TurboMole",
##    "unixyz": "UniChemXYZ",
##    "vasp": "VASP",
##    "vmol": "Viewmol",
##    "xed": "XED", # Poor description
##    "xyz": "XYZ",
##    "yob": "Yasara",
##    "zin": "Zindo",
##    "zmat": "ZMatrix"
##}
##informats = dict([(x,y) for (x,y) in _formatlookup.iteritems()
##                  if y in _informats])
##outformats = dict([(x,y) for (x,y) in _formatlookup.iteritems()
##                  if y in _outformats])

_isofact = cdk.config.IsotopeFactory.getInstance(cdk.ChemObject().getBuilder())

_bondtypes = {1: cdk.CDKConstants.BONDORDER_SINGLE,
              2: cdk.CDKConstants.BONDORDER_DOUBLE,
              3: cdk.CDKConstants.BONDORDER_TRIPLE}
_revbondtypes = dict([(y,x) for (x,y) in _bondtypes.iteritems()])

def readfile(format, filename):
    """Iterate over the molecules in a file.

    Required parameters:
       format
       filename

    You can access the first molecule in a file using:
        mol = readfile("smi", "myfile.smi").next()
        
    You can make a list of the molecules in a file using:
        mols = [mol for mol in readfile("smi", "myfile.smi")]
        
    You can iterate over the molecules in a file as shown in the
    following code snippet...

    >>> atomtotal = 0
    >>> for mol in readfile("sdf","head.sdf"):
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    if not os.path.isfile(filename):
        raise IOError, "No such file: '%s'" % filename
    builder = cdk.DefaultChemObjectBuilder.getInstance()
    if format=="sdf":
        return iter([Molecule(mol) for mol in cdk.io.iterator.IteratingMDLReader(
            java.io.FileInputStream(java.io.File(filename)),
            builder
            )])
    elif format=="smi":
        return iter([Molecule(mol) for mol in cdk.io.iterator.IteratingSmilesReader(
            java.io.FileInputStream(java.io.File(filename)),
            builder
            )])
    elif format in informats:
        reader = _informats[informats(format)]
        return iter([Molecule(reader(
            open(filename),
            builder
            ))])
    else:
        raise ValueError,"%s is not a recognised CDK format" % format

def readstring(format, string):
    """Read in a molecule from a string.

    Required parameters:
       format
       string

    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi",input)
    >>> len(mymol.atoms)
    5
    """
    if format=="smi":
        sp = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
        try:
            ans = sp.parseSmiles(string)
        except cdk.exception.InvalidSmilesException, ex:
            raise IOError, ex
        return Molecule(ans)
    elif format in informats:
        reader = _informats[informats(format)]
        return Molecule(reader(
            java.io.StringReader(string),
            cdk.DefaultChemObjectBuilder.getInstance()
            ))
    else:
        raise ValueError,"%s is not a recognised CDK format" % format

class Outputfile(object):
    """Represent a file to which *output* is to be sent.
    
    Although it's possible to write a single molecule to a file by
    calling the write() method of a molecule, if multiple molecules
    are to be written to the same file you should use the Outputfile
    class.
    
    Required parameters:
       format
       filename
    Optional parameters:
       overwrite (default is False) -- if the output file already exists,
                                       should it be overwritten?
    Methods:
       write(molecule)
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % self.filename
        if not format in outformats:
            raise ValueError,"%s is not a recognised CDK format" % format
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
        self._molwriter.write(molecule.Molecule)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.filename = None
        self._writer.close()
        self._molwriter.close()

    
class Molecule(object):
    """Represent a Pybel molecule.

    Optional parameters:
       Molecule -- a CDK Molecule (default is None)
    
    An empty Molecule is created if an Open Babel molecule is not provided.
    
    Attributes:
       atoms, charge, data, dim, energy, exactmass, flags, formula, 
       mod, molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       write(), calcfp(), calcdesc()
      
    The original CDK Molecule can be accessed using the attribute:
       Molecule
    """
    _getmethods = {
        'conformers':'GetConformers',
        # 'coords':'GetCoordinates', you can access the coordinates the atoms elsewhere
        # 'data':'GetData', has been removed
        'dim':'GetDimension',
        'energy':'GetEnergy',
        'exactmass':'GetExactMass',
        'flags':'GetFlags',
        'formula':'GetFormula',
        # 'internalcoord':'GetInternalCoord', # Causes SWIG warning
        'mod':'GetMod',
        'molwt':'GetMolWt',
        'sssr':'GetSSSR',
        'title':'GetTitle',
        'charge':'GetTotalCharge',
        'spin':'GetTotalSpinMultiplicity'
    }
    
    def __init__(self, Molecule):
        
        if hasattr(Molecule, "_xchange"):
            Molecule = readstring("smi", Molecule._xchange).Molecule
        self.Molecule = Molecule
        
    def __getattr__(self, attr):
        """Return the value of an attribute

        Note: The values are calculated on-the-fly. You may want to store the value in
        a variable if you repeatedly access the same attribute.
        """
        if attr == "atoms":
            return [Atom(self.Molecule.getAtom(i)) for i in range(self.Molecule.getAtomCount())]
        elif attr == "_atoms":
            ans = []
            for i in range(self.Molecule.getAtomCount()):
                atom = self.Molecule.getAtom(i)
                _isofact.configure(atom)
                ans.append( (atom.getAtomicNumber(),) )
            return ans
##        elif attr == 'exactmass':
              # Is supposed to use the most abundant isotope but
              # actually uses the next most abundant
##            return cdk.tools.MFAnalyser(self.Molecule).getMass()
        elif attr == "data":
            return MoleculeData(self.Molecule)
        elif attr == 'molwt':
            return cdk.tools.MFAnalyser(self.Molecule).getCanonicalMass()
        elif attr == 'formula':
            return cdk.tools.MFAnalyser(self.Molecule).getMolecularFormula()
        elif attr == "_bonds":
            ans = []
            for i in range(self.Molecule.getBondCount()):
                bond = self.Molecule.getBond(i)
                bo = bond.getOrder()
                atoms = [self.Molecule.getAtomNumber(x) for x in bond.atoms()]
                ans.append( (atoms[0], atoms[1], _revbondtypes[bo]) )
            return ans
        elif attr == "_exchange":
            return self.write("smi")
        else:
            raise AttributeError, "Molecule has no attribute '%s'" % attr

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.
        
        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        for atom in self.atoms:
            yield atom

    def __str__(self):
        return self.write()

    def addh(self):
        hAdder = cdk.tools.HydrogenAdder()
        hAdder.addExplicitHydrogensToSatisfyValency(self.Molecule)

    def removeh(self):
        atommanip = cdk.tools.manipulator.AtomContainerManipulator()
        atommanip.removeHydrogens(self.Molecule)

    def write(self, format="smi", filename=None, overwrite=False):       
        if format not in outformats:
            raise ValueError,"%s is not a recognised CDK format" % format
        if filename == None:
            if format == "smi":
                sg = cdk.smiles.SmilesGenerator()
                return sg.createSMILES(self.Molecule)
            else:
                writer = java.io.StringWriter()
        else:
            if not overwrite and os.path.isfile(filename):
                raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename            
            writer = java.io.FileWriter(java.io.File(filename))
        _outformats[format](writer).writeMolecule(self.Molecule)
        writer.close()
        if filename == None:
            return writer.toString()

    def __iter__(self):
        return iter(self.__getattr__("atoms"))

    def __str__(self):
        return self.write("smi")

    # def addh(self):
    #    ha = cdk.tools.HydrogenAdder()
    #    ha.addExplicitHydrogensToSatisfyValency(self.Molecule)        

    def calcfp(self, fp="daylight"):
        # if fp == "substructure":
        #    fingerprinter = cdk.fingerprint.SubstructureFingerprinter(
        #        cdk.fingerprint.StandardSubstructureSets.getFunctionalGroupSubstructureSet()
        #        )
        fp = fp.lower()
        if fp == "graph":
            fingerprinter = cdk.fingerprint.GraphOnlyFingerprinter()
        elif fp == "daylight":
            fingerprinter = cdk.fingerprint.Fingerprinter()
        else:
            raise ValueError, "%s is not a recognised CDK Fingerprint type" % fp
        return Fingerprint(fingerprinter.getFingerprint(self.Molecule))

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, the full list of Open Babel
        descriptors is calculated: LogP, PSA and MR.
        """
        if not descnames:
            descnames = descriptors
        ans = {}
        # Clone it to add hydrogens
        clone = self.Molecule.clone()
        Molecule(clone).addh()
        for descname in descnames:
            # Clone it to workaround CDK1.0.2 bug where
            # different descriptor calculations are not
            # independent (should be fixed for next release)
            cloneagain = clone.clone()
            try:
                desc = _descdict[descname]
            except KeyError:
                raise ValueError, "%s is not a recognised CDK descriptor type" % descname
            try:
                value = desc.calculate(cloneagain).getValue()
                if hasattr(value, "get"): # Instead of array
                    for i in range(value.length()):
                        ans[descname + ".%d" % i] = value.get(i)
                elif hasattr(value, "doubleValue"):
                    ans[descname] = value.doubleValue()
                else:
                    ans[descname] = value.intValue()
            except cdk.exception.CDKException, ex:
                # Can happen if molecule has no 3D coordinates
                pass
            except java.lang.NullPointerException, ex:
                # Happens with moment of inertia descriptor
                pass
        return ans    

    def draw(self, show=True, filename=None, update=False, web=False,
             usecoords=False):
        writetofile = filename is not None

        if not usecoords:            
            # Do the SDG
            sdg = cdk.layout.StructureDiagramGenerator()
            sdg.setMolecule(self.Molecule)
            sdg.generateExperimentalCoordinates()
            newmol = Molecule(sdg.getMolecule())
            if update:
                for atom, newatom in zip(self.atoms, newmol.atoms):
                    coords = newatom.Atom.getPoint2d()
                    atom.Atom.setPoint3d(javax.vecmath.Point3d(
                                         coords.x, coords.y, 0.0))
        else:
            newmol = self
            
        if writetofile or show:
            if writetofile:
                filedes = None
            else:
                filedes, filename = tempfile.mkstemp()
            if not web:
                # Create OASA molecule
                mol = oasa.molecule()
                for atom, newatom in zip(self._atoms, newmol.atoms):
                    if not usecoords:
                        coords = newatom.Atom.getPoint2d()
                    else:
                        coords = newatom.Atom.getPoint3d()
                    v = mol.create_vertex()
                    v.symbol = _isofact.getElement(atom[0]).getSymbol()
                    mol.add_vertex(v)
                    v.x, v.y, v.z = coords.x * 30., coords.y * -30., 0.0
                for bond in self._bonds:
                    e = mol.create_edge()
                    e.order = bond[2]
                    mol.add_edge(bond[0], bond[1], e)                        
                oasa.cairo_out.cairo_out().mol_to_cairo(mol, filename)
            else:
                encodesmiles = base64.urlsafe_b64encode(bz2.compress(self.write("smi")))
                imagedata = urllib.urlopen("http://www.chembiogrid.org/cheminfo/rest/depict/" +
                                       encodesmiles).read()
                if writetofile:
                    print >> open(filename, "wb"), imagedata
            if show:
                root = tk.Tk()
                root.title((hasattr(self, "title") and self.title)
                           or self.__str__().rstrip())
                frame = tk.Frame(root, colormap="new", visual='truecolor').pack()
                if web:
                    image = PIL.open(StringIO.StringIO(imagedata))
                else:
                    image = PIL.open(filename)
                imagedata = piltk.PhotoImage(image)
                label = tk.Label(frame, image=imagedata).pack()
                quitbutton = tk.Button(root, text="Close", command=root.destroy).pack(fill=tk.X)
                root.mainloop()
            if filedes:
                os.close(filedes)
                os.remove(filename)

##    def localopt(self, forcefield="MMFF94", steps=100):
##        forcefield = forcefield.lower()
##        
##        geoopt = cdk.modeling.forcefield.GeometricMinimizer()
##        geoopt.setMolecule(self.Molecule, False)
##        points = [t.Atom.point3d for t in self.atoms]
##        coords = [(t.x, t.y, t.z) for t in points]
##        print coords
##        if forcefield == "MMFF94":
##            ff = cdk.modeling.forcefield.MMFF94Energy(self.Molecule,
##                                geoopt.getPotentialParameterSet())
##        # geoopt.setMMFF94Tables()
##        geoopt.steepestDescentsMinimization(coords, forcefield)
##
##    def make3D(self, forcefield="MMFF94", steps=50):
##        """Generate 3D coordinates.
##        
##        Optional parameters:
##           forcefield -- default is "MMFF94"
##           steps -- default is 50
##
##        Hydrogens are added, coordinates are generated and a quick
##        local optimization is carried out with 50 steps and the
##        MMFF94 forcefield. Call localopt() if you want
##        to improve the coordinates further.
##        """
##        forcefield = forcefield.lower()
##        if forcefield not in forcefields:
##            print "hey"
##            pass
##        self.addh()
##        th3d = cdk.modeling.builder3d.TemplateHandler3D.getInstance()
##        mb3d = cdk.modeling.builder3d.ModelBuilder3D.getInstance(
##            th3d, forcefield)
##        self.Molecule = mb3d.generate3DCoordinates(self.Molecule, False)
##        self.localopt(forcefield, steps)

class Fingerprint(object):
    """A Molecular Fingerprint.
    
    Required parameters:
       obFingerprint -- a vector calculated by OBFingerprint.FindFingerprint()

    Attributes:
       fp -- the original obFingerprint
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
    """Represent a Pybel atom.

    Optional parameters:
       OBAtom -- an Open Babel Atom (default is None)
       index -- the index of the atom in the molecule (default is None)
     
    An empty Atom is created if an Open Babel atom is not provided.
    
    Attributes:
       atomicmass, atomicnum, cidx, coords, coordidx, exactmass,
       formalcharge, heavyvalence, heterovalence, hyb, idx,
       implicitvalence, index, isotope, partialcharge, spin, type,
       valence, vector.

    (refer to the Open Babel library documentation for more info).
    
    The original Open Babel atom can be accessed using the attribute:
       OBAtom
    """
    
    _getmethods = {
        'atomicmass':'GetAtomicMass',
        'atomicnum':'GetAtomicNum',
        'cidx':'GetCIdx',
        'coordidx':'GetCoordinateIdx',
        # 'data':'GetData', has been removed
        'exactmass':'GetExactMass',
        'formalcharge':'GetFormalCharge',
        'heavyvalence':'GetHvyValence',
        'heterovalence':'GetHeteroValence',
        'hyb':'GetHyb',
        'idx':'GetIdx',
        'implicitvalence':'GetImplicitValence',
        'isotope':'GetIsotope',
        'partialcharge':'GetPartialCharge',
        'spin':'GetSpinMultiplicity',
        'type':'GetType',
        'valence':'GetValence',
        'vector':'GetVector',
        }

    def __init__(self, Atom):
        self.Atom = Atom
        
    def __getattr__(self, attr):
        if attr == "coords":
            coords = self.Atom.point3d
            if not coords:
                return (0., 0., 0.)
            else:
                return (coords.x, coords.y, coords.z)
        elif attr == "atomicnum":
            _isofact.configure(self.Atom)
            return self.Atom.getAtomicNumber()
        else:
            raise AttributeError, "Atom has no attribute %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
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
        self.smarts = cdk.smiles.smarts.SMARTSQueryTool(smartspattern)
        
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
      obmol -- an Open Babel OBMol 

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying OBMol.

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
    def __init__(self, mol):
        self._mol = mol
    def _data(self):
        return self._mol.getProperties()
    def _testforkey(self, key):
        if not key in self:
            raise KeyError, "'%s'" % key
    def keys(self):
        return list(self._data().keys())
    def values(self):
        return list(self._data().values())
    def items(self):
        return [(k, self[k]) for k in self._data().keys()]
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

##>>> readstring("smi", "CCC").calcfp().bits
##[542, 637, 742]
##>>> readstring("smi", "CCC").calcfp("graph").bits
##[595, 742, 927]
##>>> readstring("smi", "CC=C").calcfp().bits
##[500, 588, 637, 742]
##>>> readstring("smi", "CC=C").calcfp("graph").bits
##[595, 742, 927]
##>>> readstring("smi", "C").calcfp().bits
##[742]
##>>> readstring("smi", "CC").calcfp().bits
##[637, 742]

# >>> readstring("smi", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# CCCCCCCCCCCCCC").calcdesc(["lipinskifailures"])
# {'lipinskifailures': 1}

if __name__=="__main__": #pragma: no cover
    mol = readstring("smi", "CCCC")
    print mol

    for mol in readfile("sdf", "head.sdf"):
        pass
    #mol = readstring("smi","CCN(CC)CC") # triethylamine
    #smarts = Smarts("[#6][#6]")
    # print smarts.findall(mol)
    mol = readstring("smi", "CC=O")
    # d = mol.calcdesc()
