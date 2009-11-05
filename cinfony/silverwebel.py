"""
webel
"""
import re
from time import sleep

# .NET classes
from System.Net import WebClient
from System import Uri, UriKind
_webclient = WebClient()

tk = None

informats = ["smi", "inchikey", "inchi", "name"]
"""A dictionary of supported input formats"""
outformats = ["smi", "cdxml", "inchi", "sdf", "names", "inchikey",
              "alc", "cerius", "charmm", "cif", "cml", "ctx",
              "gjf", "gromacs", "hyperchem", "jme", "maestro",
              "mol", "mol2", "mrv", "pdb", "sdf3000", "sln", "xyz"]
"""A dictionary of supported output formats"""

fps = ["std", "maccs", "estate"]
"""A list of supported fingerprint types"""

def _esc(text):
    return text.replace("#", "%23")
# The following function is taken from urllib.py in the IronPython dist
def _quo(text, safe="/"):
    always_safe = ('ABCDEFGHIJKLMNOPQRSTUVWXYZ'
               'abcdefghijklmnopqrstuvwxyz'
               '0123456789' '_.-')
    _safemaps = {}
    cachekey = (safe, always_safe)
    try:
        safe_map = _safemaps[cachekey]
    except KeyError:
        safe += always_safe
        safe_map = {}
        for i in range(256):
            c = chr(i)
            safe_map[c] = (c in safe) and c or ('%%%02X' % i)
        _safemaps[cachekey] = safe_map
    res = map(safe_map.__getitem__, text)
    return ''.join(res)

def _makeserver(serverurl):
    """Curry the name of the server"""
    def server(*urlcomponents):
        url = "/%s/" % serverurl + "/".join(urlcomponents)
        result = [False, None, None]
        def callback(s, e):
            result[0] = True
            result[1] = e.Error
            if not result[1]:
                result[2] = e.Result
        webclient = WebClient()
        webclient.DownloadStringCompleted += callback
        webclient.DownloadStringAsync(Uri(url, UriKind.Relative))
        while not result[0]:
            sleep(0.5)
        if result[1]:
            raise IOError, "Problem accessing web server\n%s" % result[1]
        return result[2]
    return server

rajweb = _makeserver("rajweb")
nci = _makeserver("nci")

_descs = None # Cache the list of descriptors
def getdescs():
    """Return a list of supported descriptor types"""
    global _descs
    if not _descs:
        response = rajweb("descriptors").rstrip()
        _descs = [x.split(".")[-1] for x in response.split("\n")]
    return _descs

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
    if not format in informats:
        raise ValueError("%s is not a recognised Webel format" % format)
    
    if format != "smi":
        smiles = nci(string, "smiles").rstrip()
    else:
        smiles = string
    mol = Molecule(smiles)
    if format == "name":
        mol.title = string
    return mol

class Outputfile(object):
    """Represent a file to which *output* is to be sent.
   
    Although it's possible to write a single molecule to a file by
    calling the write() method of a molecule, if multiple molecules
    are to be written to the same file you should use the Outputfile
    class.
    
    Required parameters:
       format - see the outformats variable for a list of available
                output formats
       filename

    Optional parameters:
       overwrite -- if the output file already exists, should it
                   be overwritten? (default is False)
                   
    Methods:
       write(molecule)
       close()
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format.lower()
        if not overwrite and os.path.isfile(filename):
            raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % self.filename)
        if not format in outformats:
            raise ValueError("%s is not a recognised Webel format" % format)
        self.file = open(filename, "w")
    
    def write(self, molecule):
        """Write a molecule to the output file.
        
        Required parameters:
           molecule
        """
        if self.file.closed:
            raise IOError("Outputfile instance is closed.")
        if self.format == "smi":
            output = molecule.smiles
        else:
            output = nci(_esc(molecule.smiles), "file?format=%s" % self.format).rstrip()
        
        print >> self.file, output

    def close(self):
        """Close the Outputfile to further writing."""
        self.file.close()

class Molecule(object):
    """Represent a Pybel Molecule.

    Required parameter:
       OBMol -- an Open Babel OBMol or any type of cinfony Molecule
 
    Attributes:
       atoms, charge, conformers, data, dim, energy, exactmass, formula, 
       molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       addh(), calcfp(), calcdesc(), draw(), localopt(), make3D(), removeh(),
       write() 
      
    The underlying Open Babel molecule can be accessed using the attribute:
       OBMol
    """
    _cinfony = True

    def __init__(self, smiles):
        
##        if hasattr(OBMol, "_cinfony"):
##            a, b = OBMol._exchange
##            if a == 0:
##                mol = readstring("smi", b)
##            else:
##                mol = readstring("mol", b)
##            OBMol = mol.OBMol

        self.smiles = smiles
        self.title = ""
 
    @property
    def formula(self): return rajweb("mf", _quo(self.smiles))
    @property
    def molwt(self): return float(rajweb("mw", _quo(self.smiles)))
##    @property
##    def _exchange(self):
##        if self.OBMol.HasNonZeroCoords():
##            return (1, self.write("mol"))
##        else:
##            return (0, self.write("can").split()[0])

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, all available descriptors are
        calculated. See the descs variable for a list of available
        descriptors.
        """
        if not descnames:
            descnames = getdescs()
        else:
            for descname in descnames:
                if descname not in getdescs():
                    raise ValueError("%s is not a recognised Webel descriptor type" % descname)
        ans = {}
        p = re.compile("""Descriptor parent="(\w*)" name="([\w\-\+\d]*)" value="([\d\.]*)""")
        for descname in descnames:
            longname = "org.openscience.cdk.qsar.descriptors.molecular." + descname
            response = rajweb("descriptor", longname, _quo(self.smiles))
            for match in p.findall(response):
                if match[2]:
                    ans["%s_%s" % (match[0], match[1])] = float(match[2])
        return ans
    
    def calcfp(self, fptype="std"):
        """Calculate a molecular fingerprint.
        
        Optional parameters:
           fptype -- the fingerprint type (default is "std"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        fptype = fptype.lower()
        if fptype not in fps:
            raise ValueError("%s is not a recognised Webel Fingerprint type" % fptype)
        fp = rajweb("fingerprint/%s/%s" % (fptype, _quo(self.smiles))).rstrip()
        return Fingerprint(fp)

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
        if not format in outformats:
            raise ValueError("%s is not a recognised Webel format" % format)
        if format == "smi":
            output = self.smiles
        elif format == "names":
            try:
                output = nci(_esc(self.smiles), "%s" % format).rstrip().split("\n")
            except urllib2.URLError, e:
                if e.code == 404:
                    output = []
        else:
            output = nci(_esc(self.smiles), "file?format=%s" % format).rstrip()

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % filename)
            outputfile = open(filename, "w")
            print >> outputfile, output
            outputfile.close()
        else:
            return output

##    def make3D(self, forcefield = "mmff94", steps = 50):
##        """Generate 3D coordinates.
##        
##        Optional parameters:
##           forcefield -- default is "mmff94". See the forcefields variable
##                         for a list of available forcefields.
##           steps -- default is 50
##
##        Once coordinates are generated, hydrogens are added and a quick
##        local optimization is carried out with 50 steps and the
##        MMFF94 forcefield. Call localopt() if you want
##        to improve the coordinates further.
##        """
##        forcefield = forcefield.lower()
##        _builder.Build(self.OBMol)
##        self.addh()
##        self.localopt(forcefield, steps)
##
    def __str__(self):
        return self.write()

    def draw(self, show=True, filename=None):
        """Create a 2D depiction of the molecule.

        Optional parameters:
          show -- display on screen (default is True)
          filename -- write to file (default is None)

        Tkinter and Python Imaging Library are required for
        image display.
        """
        imagedata = nci(_esc(self.smiles), "image")
        if filename:
            print >> open(filename, "wb"), imagedata
        if show:
            if not tk:
                errormessage = ("Tkinter or Python Imaging "
                                "Library not found, but is required for image "
                                "display. See installation instructions for "
                                "more information.")
                raise ImportError, errormessage                 
            root = tk.Tk()
            root.title(self.smiles)
            frame = tk.Frame(root, colormap="new", visual='truecolor').pack()
            image = PIL.open(StringIO.StringIO(imagedata))
            imagedata = piltk.PhotoImage(image)
            label = tk.Label(frame, image=imagedata).pack()
            quitbutton = tk.Button(root, text="Close", command=root.destroy).pack(fill=tk.X)
            root.mainloop()

class Fingerprint(object):
    """A Molecular Fingerprint.
    
    Required parameters:
       fingerprint -- a vector calculated by OBFingerprint.FindFingerprint()

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
        mybits = set(self.bits)
        otherbits = set(other.bits)
        return len(mybits&otherbits) / float(len(mybits|otherbits))
    @property
    def bits(self):
        return [i for i,x in enumerate(self.fp) if x=="1"]
    def __str__(self):
        return self.fp

class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern
    
    Methods:
       findall(molecule)
    
    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> smarts.match(mol) 
    True

    The numbers returned are the indices (starting from 1) of the atoms
    that match the SMARTS pattern. In this case, there are three matches
    for each of the three ethyl groups in the molecule.
    """
    def __init__(self, smartspattern):
        """Initialise with a SMARTS pattern."""
        self.pat = smartspattern
    def match(self, molecule):
        """Does a SMARTS pattern match a particular molecule?
        
        Required parameters:
           molecule
        """
        resp = rajweb("substruct", _quo(molecule.smiles), _quo(self.pat)).rstrip()
        return resp == "true"
 
if __name__=="__main__": #pragma: no cover
    import doctest
    doctest.run_docstring_examples(rajweb, globals())
