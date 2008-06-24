from __future__ import generators

import math
import os.path
import tempfile

import org.openbabel as ob

import java.lang.System
java.lang.System.loadLibrary("openbabel")

def _formatstodict(list):
    list = [list.get(i) for i in range(list.size())]
    broken = [x.replace("[Read-only]", "").replace("[Write-only]","").split(" -- ") for x in list]
    broken = [(x,y.strip()) for x,y in broken]
    return dict(broken)
_obconv = ob.OBConversion()
informats = _formatstodict(_obconv.GetSupportedInputFormat())
outformats = _formatstodict(_obconv.GetSupportedOutputFormat())

def _getplugins(findplugin, names):
    plugins = dict([(x, findplugin(x)) for x in names if findplugin(x)])
    return plugins

descs = ['LogP', 'MR', 'TPSA']
_descdict = _getplugins(ob.OBDescriptor.FindType, descs)
fps = ['FP2', 'FP3', 'FP4']
_fingerprinters = _getplugins(ob.OBFingerprint.FindFingerprint, fps)
forcefields = ['UFF', 'MMFF94', 'Ghemical']
_forcefields = _getplugins(ob.OBForceField.FindType, forcefields)
operations = ['Gen3D']
_operations = _getplugins(ob.OBOp.FindType, operations)

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
    obconversion = ob.OBConversion()
    formatok = obconversion.SetInFormat(format)
    if not formatok:
        raise ValueError,"%s is not a recognised OpenBabel format" % format
    if not os.path.isfile(filename):
        raise IOError, "No such file: '%s'" % filename
    obmol = ob.OBMol()
    notatend = obconversion.ReadFile(obmol,filename)
    while notatend:
        yield Molecule(obmol)
        obmol = ob.OBMol()
        notatend = obconversion.Read(obmol)

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
    obmol = ob.OBMol()
    obconversion = ob.OBConversion()

    formatok = obconversion.SetInFormat(format)
    if not formatok:
        raise ValueError,"%s is not a recognised OpenBabel format" % format

    success = obconversion.ReadString(obmol, string)
    if not success:
        raise IOError, "Failed to convert '%s' to format '%s'" % (
            string, format)
    return Molecule(obmol)

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
       write(molecule), close()
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % self.filename
        self.obConversion = ob.OBConversion()
        formatok = self.obConversion.SetOutFormat(self.format)
        if not formatok:
            raise ValueError,"%s is not a recognised OpenBabel format" % format
        self.total = 0 # The total number of molecules written to the file
    
    def write(self, molecule):
        """Write a molecule to the output file.
        
        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError, "Outputfile instance is closed."

        if self.total==0:
            self.obConversion.WriteFile(molecule.OBMol, self.filename)
        else:
            self.obConversion.Write(molecule.OBMol)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.obConversion.CloseOutFile()
        self.filename = None

class Molecule(object):
    """Represent a Pybel molecule.

    Required parameter:
       OBMol -- an Open Babel OBMol
       or
       Molecule -- any type of cinfony Molecule (e.g. one from cinfony.rdkit)

    If a cinfony Molecule is provided it will be converted into a pybel Molecule.       
    
    Attributes:
       atoms, charge, data, dim, energy, exactmass, flags, formula, 
       mod, molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       addh(), calcfp(), calcdesc(), draw(), localopt(), make3D(), removeh(),
       write() 
      
    The original Open Babel molecule can be accessed using the attribute:
       OBMol
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
    _cinfony = True

    def __init__(self, OBMol):
        
        if hasattr(OBMol, "_cinfony"):
            a, b = OBMol._exchange
            if a == 0:
                mol = readstring("smi", b)
            else:
                mol = readstring("mol", b)
            OBMol = mol.OBMol

        self.OBMol = OBMol
 
    def __getattr__(self, attr):
        """Return the value of an attribute

        Note: The values are calculated on-the-fly. You may want to store the value in
        a variable if you repeatedly access the same attribute.
        """
        # This function is not accessed in the case of OBMol
        if attr == "atoms":
            # Create an atoms attribute on-the-fly
            return [ Atom(self.OBMol.GetAtom(i+1)) for i in range(self.OBMol.NumAtoms()) ]
        elif attr == "data":
            # Create a data attribute on-the-fly
            return MoleculeData(self.OBMol)
        elif attr == "unitcell":
            # Create a unitcell attribute on-the-fly
            unitcell = self.OBMol.GetData(ob.openbabelConstants.UnitCell)
            if unitcell:
                return ob.openbabel.toUnitCell(unitcell)
            else:
                raise AttributeError, "Molecule has no attribute 'unitcell'"
        elif attr in self._getmethods:
            # Call the OB Method to find the attribute value
            return getattr(self.OBMol, self._getmethods[attr])()
        elif attr == "_exchange":
            if self.OBMol.HasNonZeroCoords():
                return (1, self.write("mol"))
            else:
                return (0, self.write("can").split()[0])
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

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, the full list of Open Babel
        descriptors is calculated. See pybel.descs for a list
        of available descriptors.
        """
        if not descnames:
            descnames = descs
        ans = {}
        for descname in descnames:
            try:
                desc = _descdict[descname]
            except KeyError:
                raise ValueError, "%s is not a recognised Open Babel descriptor type" % descname
            ans[descname] = desc.Predict(self.OBMol)
        return ans
    
    def calcfp(self, fptype="FP2"):
        """Calculate a molecular fingerprint.
        
        Optional parameters:
           fptype -- the name of the Open Babel fingerprint type.

        If fptype is not specified, the FP2 fingerprint type
        is used. See the Open Babel library documentation for more
        details.
        """
        fp = ob.vectorUnsignedInt()
        try:
            fingerprinter = _fingerprinters[fptype]
        except KeyError:
            raise ValueError, "%s is not a recognised Open Babel Fingerprint type" % fptype
        fingerprinter.GetFingerprint(self.OBMol, fp)
        return Fingerprint(fp)

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
        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat(format)
        if not formatok:
            raise ValueError,"%s is not a recognised OpenBabel format" % format

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename
            obconversion.WriteFile(self.OBMol,filename)
            obconversion.CloseOutFile()
        else:
            return obconversion.WriteString(self.OBMol)

    def localopt(self, forcefield="MMFF94", steps=500):
        """Locally optimize the coordinates.
        
        Optional parameters:
           forcefield -- default is "MMFF94"
           steps -- default is 500

        If the molecule does not have 2D or 3D coordinates, make3D() is
        called before the optimization. Note that the molecule needs
        to have explicit hydrogens. If not, call addh().
        """
        
        if not (self.OBMol.Has2D() or self.OBMol.Has3D()):
            self.make3D()
        ff = _forcefields[forcefield]
        ff.Setup(self.OBMol)
        ff.SteepestDescent(steps)
        ff.GetCoordinates(self.OBMol)
    
##    def globalopt(self, forcefield="MMFF94", steps=1000):
##        if not (self.OBMol.Has2D() or self.OBMol.Has3D()):
##            self.make3D()
##        self.localopt(forcefield, 250)
##        ff = _forcefields[forcefield]
##        numrots = self.OBMol.NumRotors()
##        if numrots > 0:
##            ff.WeightedRotorSearch(numrots, int(math.log(numrots + 1) * steps))
##        ff.GetCoordinates(self.OBMol)
    
    def make3D(self, forcefield="MMFF94", steps=50):
        """Generate 3D coordinates.
        
        Optional parameters:
           forcefield -- default is "MMFF94"
           steps -- default is 50

        Once coordinates are generated, hydrogens are added and a quick
        local optimization is carried out with 50 steps and the
        MMFF94 forcefield. Call localopt() if you want
        to improve the coordinates further.
        """
        _operations['Gen3D'].Do(self.OBMol)
        self.addh()
        self.localopt(forcefield, steps)

    def addh(self):
        """Add hydrogens."""
        self.OBMol.AddHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.OBMol.DeleteHydrogens()
        
    def __str__(self):
        return self.write()

class Atom(object):
    """Represent a Pybel atom.

    Required parameter:
       OBAtom -- an Open Babel OBAtom
        
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

    def __init__(self, OBAtom):
        self.OBAtom = OBAtom
        
    def __getattr__(self, attr):
        if attr == "coords":
            return (self.OBAtom.GetX(), self.OBAtom.GetY(), self.OBAtom.GetZ())
        elif attr in self._getmethods:
            return getattr(self.OBAtom, self._getmethods[attr])()
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

def _findbits(fp, bitsperint):
    """Find which bits are set in a list/vector.

    This function is used by the Fingerprint class.

    >>> _findbits([13, 71], 8)
    [1, 3, 4, 9, 10, 11, 15]
    """
    ans = []
    start = 1
    fp = [fp.get(i) for i in range(fp.size())]
    for x in fp:
        i = start
        while x > 0:
            if x % 2:
                ans.append(i)
            x >>= 1
            i += 1
        start += bitsperint
    return ans
        
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
        return ob.OBFingerprint.Tanimoto(self.fp, other.fp)
    def __getattr__(self, attr):
        if attr == "bits":
            # Create a bits attribute on-the-fly
            return _findbits(self.fp, ob.OBFingerprint.Getbitsperint())
        else:
            raise AttributeError, "Fingerprint has no attribute %s" % attr
    def __str__(self):
        return ", ".join([str(self.fp.get(i)) for i in range(self.fp.size())])

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
    [(1, 2), (4, 5), (6, 7)]

    The numbers returned are the indices (starting from 1) of the atoms
    that match the SMARTS pattern. In this case, there are three matches
    for each of the three ethyl groups in the molecule.
    """
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        self.obsmarts = ob.OBSmartsPattern()
        success = self.obsmarts.Init(smartspattern)
        if not success:
            raise IOError, "Invalid SMARTS pattern"
    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.
        
        Required parameters:
           molecule
        """
        self.obsmarts.Match(molecule.OBMol)
        vector = self.obsmarts.GetUMapList()
        return [vector.get(i) for i in range(vector.size())]
        
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
    def __init__(self, obmol):
        self._mol = obmol
    def _data(self):
        data = self._mol.GetData()
        data = [data.get(i) for i in range(data.size())]
        return [ob.openbabel.toPairData(x) for x in data
                if x.GetDataType()==ob.openbabelConstants.PairData or
                x.GetDataType()==ob.openbabelConstants.CommentData]
    def _testforkey(self, key):
        if not key in self:
            raise KeyError, "'%s'" % key
    def keys(self):
        return [x.GetAttribute() for x in self._data()]
    def values(self):
        return [x.GetValue() for x in self._data()]
    def items(self):
        return zip(self.keys(), self.values())
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self._data())
    def __contains__(self, key):
        return self._mol.HasData(key)
    def __delitem__(self, key):
        self._testforkey(key)
        self._mol.DeleteData(self._mol.GetData(key))
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
        return ob.openbabel.toPairData(self._mol.GetData(key)).GetValue()
    def __setitem__(self, key, value):
        if key in self:
            pairdata = ob.openbabel.toPairData(self._mol.GetData(key))
            pairdata.SetValue(str(value))
        else:
            pairdata = ob.OBPairData()
            pairdata.SetAttribute(key)
            pairdata.SetValue(str(value))
            # For the following line to work, you need
            #    python.security.respectJavaAccessibility = false
            # in the Jython registry
            pairdata.swigCMemOwn = 0 # So that SWIG Proxy will not delete pairdata
            self._mol.SetData(pairdata)
    def __repr__(self):
        return dict(self.iteritems()).__repr__()
 
if __name__=="__main__": #pragma: no cover
    import doctest
    doctest.testmod(verbose=True)
