import os

from Chem import AllChem
from Chem.Draw import MolDrawing
from Chem.AvailDescriptors import descDict

import DataStructs
import Chem.MACCSkeys
import Chem.AtomPairs.Pairs
import Chem.AtomPairs.Torsions

# PIL and Tkinter
try:
    import Tkinter as tk
    import Image as PIL
    import ImageTk as PILtk
except:
    PILtk = None

# Aggdraw
try:
    import aggdraw
except ImportError:
    aggdraw = None

fps = ['Daylight', 'MACCS', 'atompairs', 'torsions']
descs = descDict.keys()

_formats = {'smi': "SMILES", 'iso': "Isomeric SMILES",
            'mol': "MDL MOL file", 'sdf': "MDL SDF file"}
informats = dict([(x, _formats[x]) for x in ['mol', 'sdf', 'smi']])
outformats = dict([(x, _formats[x]) for x in ['mol', 'sdf', 'smi', 'iso']])

_forcefields = {'uff': AllChem.UFFOptimizeMolecule}
forcefields = _forcefields.keys()

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
        iterator = Chem.SDMolSupplier(filename)
        for mol in iterator:
            yield Molecule(mol)
    elif format=="mol":
        yield Molecule(Chem.MolFromMolFile(filename))
    else:
        raise ValueError,"%s is not a recognised RDKit format" % format

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
        return Molecule(Chem.MolFromMolBlock(string))
    elif format=="smi":
        mol = Chem.MolFromSmiles(string)        
    else:
        raise ValueError,"%s is not a recognised RDKit format" % format
    if mol:
        return Molecule(mol)
    else:
        raise IOError, "Failed to convert '%s' to format '%s'" % (
            string, format)

class Molecule(object):
    """Represent an RDKit molecule.

    Required parameter:
       Mol -- an RDKit Mol
       or
       Molecule -- any type of cinfony Molecule (e.g. one from cinfony.pybel)

    If a cinfony Molecule is provided it will be converted into
    an rdkit Molecule.       
    
    Attributes:
       atoms, data, molwt, title
    
    Methods:
       addh(), calcfp(), calcdesc(), draw(), localopt(), make3D(), removeh(),
       write() 
      
    The underlying RDKit Mol can be accessed using the attribute:
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
            return descDict['MolWt'](self.Mol)
        elif attr == "_exchange":
            if self.Mol.GetNumConformers() == 0:
                return (0, self.write("iso"))
            else:
                return (1, self.write("mol"))
        elif attr == "title":
            return self.Mol.GetProp("_Name")
        else:
            raise AttributeError, "Molecule has no attribute '%s'" % attr

    def addh(self):
        """Add hydrogens."""
        self.Mol = Chem.AddHs(self.Mol)
        
    def removeh(self):
        """Remove hydrogens."""
        self.Mol = Chem.RemoveHs(self.Mol)
        
    def write(self, format="smi", filename=None, overwrite=False):
        """Write the molecule to a file or return a string.
        
        Optional parameters:
           format -- default is "smi"
           filename -- default is None
           overwite -- default is False

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.
        The overwrite flag is ignored if a filename is not specified.
        It controls whether to allow an existing file to be overwritten.
        """
        format = format.lower()
        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename
        if format=="smi":
            result = Chem.MolToSmiles(self.Mol)
        elif format=="iso":
            result = Chem.MolToSmiles(self.Mol, isomericSmiles=True)
        elif format=="mol":
            result = Chem.MolToMolBlock(self.Mol)
        else:
            raise ValueError,"%s is not a recognised RDKit format" % format
        if filename:
            print >> open(filename, "w"), result
        else:
            return result

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

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, the full list of RDKit
        descriptors is calculated. See rdkit.descs for a list
        of available descriptors.
        """
        if not descnames:
            descnames = descs
        ans = {}
        for descname in descnames:
            try:
                desc = descDict[descname]
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
            fp = Fingerprint(Chem.DaylightFingerprint(self.Mol))
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

        Aggdraw is used for 2D depiction. Tkinter and
        Python Imaging Library are required for image display.
        """
        if usecoords:
            confId = 0
        else:
            if update:
                AllChem.Compute2DCoords(self.Mol)
                confId = 0
            else:
                confId = self.Mol.GetNumConformers()
                AllChem.Compute2DCoords(self.Mol, clearConfs = False)
        if show or filename:
            if not aggdraw:
                errormessage = ("Aggdraw not found, but is required for 2D"
                                "structure depiction. "
                                "See installation instructions for more "
                                "information.")
                raise ImportError, errormessage
              
            Chem.Kekulize(self.Mol)
            MolDrawing.registerCanvas('agg')
            img = PIL.new("RGBA",(300,300),"white")
            canvas = aggdraw.Draw(img)                
            canvas.setantialias(True)
            drawer = MolDrawing.MolDrawing(canvas)
            drawer.wedgeDashedBonds = True
            drawer.AddMol(self.Mol, confId = confId)
            canvas.flush()
        
            if filename: # Note: overwrite is allowed          
                img.save(filename)
            if show:
                if not tk:
                    errormessage = ("Tkinter or Python Imaging "
                                    "Library not found, but is required for image "
                                    "display. See installation instructions for "
                                    "more information.")
                    raise ImportError, errormessage                  
                root = tk.Tk()
                root.title((hasattr(self, "title") and self.title)
                           or self.__str__().rstrip())
                frame = tk.Frame(root, colormap="new", visual='truecolor').pack()
                imagedata = PILtk.PhotoImage(img)
                label = tk.Label(frame, image=imagedata).pack()
                quitbutton = tk.Button(root, text="Close", command=root.destroy).pack(fill=tk.X)
                root.mainloop()
            Chem.SanitizeMol(self.Mol)

    def localopt(self, forcefield = "uff", steps = 50):
        forcefield = forcefield.lower()
        _forcefields[forcefield](self.Mol, maxIters = steps)

    def make3D(self):
        success = AllChem.EmbedMolecule(self.Mol)
        if success == -1: # Failed
            success = AllChem.EmbedMolecule(self.Mol,
                                            useRandomCoords = True)
            if success == -1:
                raise Error, "Embedding failed!"
        self.localopt()
        
class Atom(object):
    """Represent a Pybel atom.

    Required parameters:
       Atom -- an RDKit Atom
     
    Attributes:
       atomicmass, atomicnum, cidx, coords, coordidx, exactmass,
       formalcharge, heavyvalence, heterovalence, hyb, idx,
       implicitvalence, index, isotope, partialcharge, spin, type,
       valence, vector.

    (refer to the Open Babel library documentation for more info).
    
    The original RDKit Atom can be accessed using the attribute:
       Atom
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

    def __init__(self, Atom):
        self.Atom = Atom
        
    def __getattr__(self, attr):
        if attr == "coords":
            owningmol = self.Atom.GetOwningMol()
            if owningmol.GetNumConformers() == 0:
                raise AttributeError, "Atom has no coordinates (0D structure)"
            idx = self.Atom.GetIdx()
            atomcoords = owningmol.GetConformer().GetAtomPosition(idx)
            return (atomcoords[0], atomcoords[1], atomcoords[2])
        elif attr in self._getmethods:
            return getattr(self.Atom, self._getmethods[attr])()        
        else:
            raise AttributeError, "Atom has no attribute %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
        if hasattr(self, "coords"):
            return "Atom: %d (%.2f %.2f %.2f)" % (self.atomicnum, self.coords[0],
                                                    self.coords[1], self.coords[2])
        else:
            return "Atom: %d (no coords)" % (self.atomicnum)

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
            self._writer = Chem.SDWriter(self.filename)
        elif format=="smi":
            self._writer = Chem.SmilesWriter(self.filename)
        else:
            raise ValueError,"%s is not a recognised RDKit format" % format
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
       findall(molecule)
    
    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol) 
    [(0, 1), (3, 4), (5, 6)]
    """
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        self.rdksmarts = Chem.MolFromSmarts(smartspattern)
        if not self.rdksmarts:
            raise IOError, "Invalid SMARTS pattern."

    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.
        
        Required parameters:
           molecule
        """
        return molecule.Mol.GetSubstructMatches(self.rdksmarts)

class MoleculeData(object):
    """Store molecule data in a dictionary-type object
    
    Required parameters:
      obmol -- an RDKit Mol 

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
    
    Required parameters: # FIXME
       fingerprint -- a vector calculated by OBFingerprint.FindFingerprint()

    Attributes:
       fp -- the original fingerprint
       bits -- a list of bits set in the Fingerprint

    Methods:
       The "|" operator can be used to calculate the Tanimoto coeff. For example,
       given two Fingerprints 'a', and 'b', the Tanimoto coefficient is given by:
          tanimoto = a | b
    """
    def __init__(self, fingerprint):
        self.fp = fingerprint
    def __or__(self, other):
        return DataStructs.FingerprintSimilarity(self.fp, other.fp)
    def __getattr__(self, attr):
        if attr == "bits":
            # Create a bits attribute on-the-fly
            return list(self.fp.GetOnBits())
        else:
            raise AttributeError, "Fingerprint has no attribute %s" % attr
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
    
