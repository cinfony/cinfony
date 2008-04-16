import os
from jpype import *

jvm = os.environ['JPYPE_JVM']
cp = os.environ['CLASSPATH']
startJVM(jvm, "-Djava.class.path=" + cp)

cdk = JPackage("org").openscience.cdk

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
##    "dock": "Dock 5 Box",
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
        raise ValueError,"%s is not a recognised OpenBabel format" % format

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
        return Molecule(sp.parseSmiles(string))
    elif format in informats:
        reader = _informats[informats(format)]
        return Molecule(reader(
            java.io.StringReader(string),
            cdk.DefaultChemObjectBuilder.getInstance()
            ))
    else:
        raise ValueError,"%s is not a recognised OpenBabel format" % format
    
class Molecule(object):
    """Represent a Pybel molecule.

    Optional parameters:
       CDKMol -- an Open Babel molecule (default is None)
    
    An empty Molecule is created if an Open Babel molecule is not provided.
    
    Attributes:
       atoms, charge, data, dim, energy, exactmass, flags, formula, 
       mod, molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       write(), calcfp(), calcdesc()
      
    The original Open Babel molecule can be accessed using the attribute:
       CDKMol
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
    
    def __init__(self, CDKMol):

        self.CDKMol = CDKMol
        
    def __getattr__(self, attr):
        """Return the value of an attribute

        Note: The values are calculated on-the-fly. You may want to store the value in
        a variable if you repeatedly access the same attribute.
        """
        if attr == "atoms":
            return [Atom(self.CDKMol.getAtom(i)) for i in range(self.CDKMol.getAtomCount()) ]
        elif attr == 'exactmass':
            # I have probably confused exact, canonical and natural masses
            return cdk.formula.MolecularFormulaManipulator.getMolecularFormula(self.CDKMol).getCanonicalMass()
        elif attr == 'molwt':
            # I have probably confused exact, canonical and natural masses
            return cdk.formula.MolecularFormulaManipulator.getTotalExactMass(
                cdk.formula.MolecularFormulaManipulator.getMolecularFormula(
                    self.CDKMol))
        elif attr == 'formula':
            return cdk.formula.MolecularFormulaManipulator.getString(
                cdk.formula.MolecularFormulaManipulator.getMolecularFormula(
                    self.CDKMol))
        else:
            raise AttributeError, "Molecule has no attribute '%s'" % attr

    def write(self, format, filename=None):
        if filename==None:
            if format=="smi":
                sg = cdk.smiles.SmilesGenerator()
                return sg.createSMILES(self.CDKMol)
        else:
            cdk.io.SMILESWriter(open(filename, "w")).writeMolecule(self.CDKMol)

    def __iter__(self):
        return iter(self.__getattr__("atoms"))

    def __str__(self):
        return self.write("smi")

    # def addh(self):
    #    ha = cdk.tools.HydrogenAdder()
    #    ha.addExplicitHydrogensToSatisfyValency(self.CDKMol)        

    def calcfp(self, fp="daylight"):
        if fp == "substructure":
            fingerprint = cdk.fingerprint.SubstructureFingerprinter
        else:
            fingerprinter = cdk.fingerprint.Fingerprinter
        return fingerprinter.getFingerprint(self.CDKMol)

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
        elif attr in self._getmethods:
            return getattr(self.OBAtom, self._getmethods[attr])()
        else:
            raise AttributeError, "Molecule has no attribute %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
        return "Atom: %d %s" % (self.atomicnum, self.coords.__str__())

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
        match = self.smarts.matches(molecule.CDKMol)
        return list(self.smarts.getUniqueMatchingAtoms())

if __name__=="__main__": #pragma: no cover
    mol = readstring("smi", "CCCC")
    print mol

    for mol in readfile("sdf", "head.sdf"):
        print mol.formula, mol.write("smi")

    mol = readstring("smi","CCN(CC)CC") # triethylamine
    smarts = Smarts("[#6][#6]")
    print smarts.findall(mol)


    
    

    
