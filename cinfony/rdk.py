import os
import cinfony

from Chem import AllChem as rdkit
from Chem.AvailDescriptors import descDict as descriptors
import Chem.MACCSkeys
import Chem.AtomPairs.Pairs
import Chem.AtomPairs.Torsions
import DataStructs

fps = ['Daylight', 'MACCS', 'atompairs', 'torsions']

_formats = {'smi': "SMILES", 'iso': "Isomeric SMILES",
            'mol': "MDL MOL file", 'sdf': "MDL SDF file"}
informats = dict([(x, _formats[x]) for x in ['mol', 'sdf', 'smi']])
outformats = dict([(x, _formats[x]) for x in ['mol', 'sdf', 'smi', 'iso']])

_bondtypes = {1: rdkit.BondType.SINGLE,
              2: rdkit.BondType.DOUBLE,
              3: rdkit.BondType.TRIPLE}
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
    ...     mol.addh()
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    if not os.path.isfile(filename):
        raise IOError, "No such file: '%s'" % filename    
    format = format.lower()    
    if format=="sdf":
        iterator = rdkit.SDMolSupplier(filename)
        for mol in iterator:
            yield Molecule(mol)
    elif format=="mol":
        yield Molecule(rdkit.MolFromMolFile(filename))
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
    format = format.lower()    
    if format=="mol":
        return Molecule(rdkit.MolFromMolBlock(string))
    elif format=="smi":
        return Molecule(rdkit.MolFromSmiles(string))
    else:
        raise ValueError,"%s is not a recognised OpenBabel format" % format
       

class Molecule(cinfony.Molecule):
    """Represent an RDKit molecule.

    Required parameter:
       Mol -- an RDKit Mol or any cinfony Molecule
    
    Attributes:
       atoms, charge, data, dim, energy, exactmass, flags, formula, 
       mod, molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       write(), calcfp(), calcdesc()
      
    The underlying RDKit Mol can be accessed using the attribute:
       Mol
    """

    def __init__(self, Mol):
        if hasattr(Mol, "atoms") and hasattr(Mol, "_bonds"):
            Mol = self._buildMol(Mol)
        self.Mol = Mol
   
    def _buildMol(self, molecule):
        rdmol = rdkit.Mol()
        rdedmol = rdkit.EditableMol(rdmol)
        for atom in molecule.atoms:
            rdatom = rdkit.Atom(atom.atomicnum)
            if hasattr(atom, "formalcharge"):
                rdatom.SetFormalCharge(atom.formalcharge)
            rdedmol.AddAtom(rdatom)
        
        for bond in molecule._bonds:
            rdedmol.AddBond(bond[0],
                            bond[1],
                            _bondtypes[bond[2]])
        rdmol = rdedmol.GetMol()
        rdkit.SanitizeMol(rdmol)
        return rdmol        

    def __getattr__(self, attr):
        """Return the value of an attribute

        Note: The values are calculated on-the-fly. You may want to store the value in
        a variable if you repeatedly access the same attribute.
        """
        # This function is not accessed in the case of OBMol
        if attr == "atoms":
            return [Atom(rdkatom) for rdkatom in self.Mol.GetAtoms()]
        elif attr == "data":
            return MoleculeData(self.Mol)
        elif attr == "molwt":
            return descriptors['MolWt'](self.Mol)
        elif attr == "_bonds":
            rdkit.Kekulize(self.Mol)
            ans = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx(),
                    _revbondtypes[x.GetBondType()])
                   for x in self.Mol.GetBonds()]
            rdkit.SanitizeMol(self.Mol)
            return ans
        else:
            raise AttributeError, "Molecule has no attribute '%s'" % attr

    def addh(self):
        self.Mol = rdkit.AddHs(self.Mol)

    def write(self, format="SMI", filename=None, overwrite=False):
        """Write the molecule to a file or return a string.
        
        Optional parameters:
           format -- default is "SMI"
           filename -- default is None
           overwite -- default is False

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.
        The overwrite flag is ignored if a filename is not specified.
        It controls whether to overwrite an existing file.
        """
        
        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename
        if format=="smi":
            result = rdkit.MolToSmiles(self.Mol)
        elif format=="iso":
            result = rdkit.MolToSmiles(self.Mol, isomericSmiles=True)
        elif format=="mol":
            result = rdk.MolToMolBlock(self.Mol)
        else:
            raise ValueError,"%s is not a recognised OpenBabel format" % format
        if filename:
            print >> open(filename, "w"), result
        else:
            return result

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, the full list of Open Babel
        descriptors is calculated: LogP, PSA and MR.
        """
        if not descnames:
            descnames = descriptors.keys()
        ans = {}
        for descname in descnames:
            try:
                desc = descriptors[descname]
            except KeyError:
                raise ValueError, "%s is not a recognised RDKit descriptor type" % descname
            ans[descname] = desc(self.Mol)
        return ans

    def calcfp(self, fptype="Daylight"):
        """Calculate a molecular fingerprint.
        
        Optional parameters:
           fptype -- the name of the Open Babel fingerprint type.

        If fptype is not specified, the FP2 fingerprint type
        is used. See the Open Babel library documentation for more
        details.
        """
        fptype = fptype.lower()
        if fptype=="daylight":
            fp = Fingerprint(rdkit.DaylightFingerprint(self.Mol))
        elif fptype=="maccs":
            fp = Fingerprint(Chem.MACCSkeys.GenMACCSKeys(self.Mol))
        elif fptype=="atompairs":
            # Going to leave as-is. See Atom Pairs documentation.
            fp = Chem.AtomPairs.Pairs.GetAtomPairFingerprintAsIntVect(self.Mol)
        elif fptype=="torsions":
            # Going to leave as-is.
            fp = Chem.AtomPairs.Torsions.GetTopologicalTorsionFingerprintAsIntVect(self.Mol)
        else:
            raise ValueError, "%s is not a recognised RDKit Fingerprint type" % fptype
        return fp

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
        # 'exactmass':'GetMass',
        'formalcharge':'GetFormalCharge',
        'heavyvalence':'GetHvyValence',
        'heterovalence':'GetHeteroValence',
        #'hyb':'GetHyb',
        #'idx':'GetIdx',
        #'implicitvalence':'GetImplicitValence',
        #'isotope':'GetIsotope',
        #'partialcharge':'GetPartialCharge',
        #'spin':'GetSpinMultiplicity',
        #'type':'GetType',
        #'valence':'GetValence',
        #'vector':'GetVector',
        }

    def __init__(self, RDKatom):
        self.RDKatom = RDKatom
        
    def __getattr__(self, attr):
        if attr == "coords":
            idx = self.RDKatom.GetIdx()
            conformer = self.RDKatom.GetOwningMol().GetConformer()
            atomcoords = conformer.GetAtomPosition(idx)
            return (atomcoords[0], atomcoords[1], atomcoords[2])
        elif attr in self._getmethods:
            return getattr(self.RDKatom, self._getmethods[attr])()        
        else:
            raise AttributeError, "Atom has no attribute %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
        return "Atom: %d (%.4f, %.4f, %.4f)" % (self.atomicnum, self.coords[0],
                                                self.coords[1], self.coords[2])


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
        if format=="sdf":
            self._writer = rdkit.SDWriter(self.filename)
        elif format=="smi":
            self._writer = rdkit.SmilesWriter(self.filename)
        else:
            raise ValueError,"%s is not a recognised OpenBabel format" % format
        self.total = 0 # The total number of molecules written to the file
    
    def write(self, molecule):
        """Write a molecule to the output file.
        
        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError, "Outputfile instance is closed."

        self._writer.write(molecule.Mol)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.filename = None
        self._writer.flush()
        del self._writer

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
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        self.rdksmarts = rdkit.MolFromSmarts(smartspattern)

    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.
        
        Required parameters:
           molecule
        """
        return molecule.Mol.GetSubstructMatches(self.rdksmarts)

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
    def __init__(self, Mol):
        self._mol = Mol
    def _testforkey(self, key):
        if not key in self:
            raise KeyError, "'%s'" % key
    def keys(self):
        return self._mol.GetPropNames()
    def values(self):
        return [self._mol.GetProp(x) for x in self.keys()]
    def items(self):
        return zip(self.keys(), self.values())
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self.keys())
    def __contains__(self, key):
        return self._mol.HasProp(key)
    def __delitem__(self, key):
        self._testforkey(key)
        self._mol.ClearProp(key)
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
        return self._mol.GetProp(key)
    def __setitem__(self, key, value):
        self._mol.SetProp(key, str(value))
    def __repr__(self):
        return dict(self.iteritems()).__repr__()

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
    def __init__(self, obFingerprint):
        self.fp = obFingerprint
    def __or__(self, other):
        return DataStructs.FingerprintSimilarity(self.fp, other.fp)
    def __getattr__(self, attr):
        if attr == "bits":
            # Create a bits attribute on-the-fly
            return list(self.fp.GetOnBits())
        else:
            raise AttributeError, "Molecule has no attribute %s" % attr
    def __str__(self):
        return ", ".join([str(x) for x in _compressbits(self.fp)])

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
    
