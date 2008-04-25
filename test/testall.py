import os
import unittest

from cinfony import pybel, cdk, rdkit

# For compatability with Python2.3
try:
    from sets import Set as set
except ImportError:
    pass

if os.path.isfile("testoutput.txt"):
    os.remove("testoutput.txt")

class TestToolkit(unittest.TestCase):
    
    def setUp(self):
        self.mols = [self.toolkit.readstring("smi", "CCCC"),
                     self.toolkit.readstring("smi", "CCCN")]
        self.head = list(self.toolkit.readfile("sdf", "head.sdf"))
        self.atom = self.head[0].atoms[1]

    def FPaccesstest(self):
        # Should raise AttributeError
        return self.mols[0].calcfp().nosuchname

    def testFPTanimoto(self):
        """Test the calculation of the Tanimoto coefficient"""
        fps = [x.calcfp() for x in self.mols]
        self.assertEqual(fps[0] | fps[1], self.tanimotoresult)
        
    def testFPstringrepr(self):
        """Test the string representation and corner cases."""
        self.assertRaises(ValueError, self.mols[0].calcfp, "Nosuchname")
        self.assertRaises(AttributeError, self.FPaccesstest)
        self.assertEqual(str(self.mols[0].calcfp()),
                         '0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1073741824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0')

    def testFPbits(self):
        """Test whether the bits are set correctly."""
        bits = [x.calcfp().bits for x in self.mols]
        self.assertEqual(len(bits[0]), self.Nbits)
        bits = [set(x) for x in bits]
        # Calculate the Tanimoto coefficient the old-fashioned way
        tanimoto = len(bits[0] & bits[1]) / float(len(bits[0] | bits[1]))
        self.assertEqual(tanimoto, self.tanimotoresult)

    def RSaccesstest(self):
        # Should raise AttributeError
        return self.mols[0].nosuchname

    def testRSformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, self.toolkit.readstring, "noel", "jkjk")
    
    def testRSgetprops(self):
        """Get the values of the properties."""
        # self.assertAlmostEqual(self.mols[0].exactmass, 58.078, 3)
        # Only OpenBabel has a working exactmass
        self.assertAlmostEqual(self.mols[0].molwt, 58.12, 2)
        self.assertEqual(len(self.mols[0].atoms), 4)
        self.assertRaises(AttributeError, self.RSaccesstest)

    def testRSconversiontoMOL(self):
        """Convert to mol"""
        as_mol = self.mols[0].write("mol")
        test = """
 OpenBabel04220815032D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0
  1  2  1  0  0  0
  2  3  1  0  0  0
  3  4  1  0  0  0
M  END
"""
        for a,b in zip(test.split("\n")[2:], as_mol.split("\n")[2:]):
            self.assertEqual(a,b)

    def testRSstringrepr(self):
        """Test the string representation of a molecule"""
        self.assertEqual(str(self.mols[0]).strip(), "CCCC")

    def testRFread(self):
        """Is the right number of molecules read from the file?"""
        self.assertEqual(len(self.mols), 2)

    def RFreaderror(self):
        mol = self.toolkit.readfile("sdf", "nosuchfile.sdf").next()

    def testRFmissingfile(self):
        """Test that reading from a non-existent file raises an error."""
        self.assertRaises(IOError, self.RFreaderror)

    def RFformaterror(self):
        mol = self.toolkit.readfile("noel", "head.sdf").next()
    
    def testRFformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, self.RFformaterror)

    def RFunitcellerror(self):
        unitcell = self.mols[0].unitcell
    
    def testRFunitcellerror(self):
        """Test that accessing the unitcell raises an error"""
        self.assertRaises(AttributeError, self.RFunitcellerror)

    def testRFconversion(self):
        """Convert to smiles"""
        as_smi = [mol.write("smi").split("\t")[0] for mol in self.mols]
        test = ['CCCC', 'CCCN']
        self.assertEqual(as_smi, test)

    def testRFsingletofile(self):
        """Test the molecule.write() method"""
        mol = self.mols[0]
        mol.write("smi", "testoutput.txt")
        test = 'CCCC'
        filecontents = open("testoutput.txt", "r").readlines()[0].split("\t")[0]
        self.assertEqual(filecontents, test)
        self.assertRaises(IOError, mol.write, "smi", "testoutput.txt")
        os.remove("testoutput.txt")
        self.assertRaises(ValueError, mol.write, "noel", "testoutput.txt")
    
    def testRFoutputfile(self):
        """Test the Outputfile class"""
        self.assertRaises(ValueError, self.toolkit.Outputfile, "noel", "testoutput.txt")
        outputfile = self.toolkit.Outputfile("sdf", "testoutput.txt")
        for mol in self.head:
            outputfile.write(mol)
        outputfile.close()
        self.assertRaises(IOError, outputfile.write, mol)
        self.assertRaises(IOError, self.toolkit.Outputfile, "sdf", "testoutput.txt")
        numdollar = len([x for x in open("testoutput.txt").readlines()
                         if x.rstrip() == "$$$$"])
        os.remove("testoutput.txt")
        self.assertEqual(numdollar, 2)

    def RFdesctest(self):
        # Should raise ValueError
        self.mols[0].calcdesc("BadDescName")

    def notestRFdesc(self):
        """Test the descriptors"""
        desc = self.mols[1].calcdesc()
        self.assertEqual(len(desc), self.Ndescs)
        self.assertAlmostEqual(desc[self.tpsaname], 26.02, 2)
        self.assertRaises(ValueError, self.RFdesctest)

    def MDaccesstest(self):
        # Should raise KeyError
        return self.head[0].data['noel']

    def testMDaccess(self):
        """Change the value of a field"""
        data = self.head[0].data
        self.assertRaises(KeyError, self.MDaccesstest)
        data['noel'] = 'testvalue'
        self.assertEqual(data['noel'], 'testvalue')
        newvalues = {'hey':'there', 'yo':1}
        data.update(newvalues)
        self.assertEqual(data['yo'], '1')
        self.assertTrue('there' in data.values())

    def testMDglobalaccess(self):
        """Check out the keys"""
        data = self.head[0].data
        # self.assert_(data.has_key('Comment'))
        self.assertFalse(data.has_key('Noel'))
        self.assertEqual(len(data), 2)
        for key in data:
            self.assertEqual(key in ['Comment', 'NSC'], True)
        self.assertEqual(repr(data), "{'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}")

    def testMDdelete(self):
        """Delete some keys"""
        data = self.head[0].data
        self.assertTrue(data.has_key('NSC'))
        del data['NSC']
        self.assertFalse(data.has_key('NSC'))
        data.clear()
        self.assertEqual(len(data), 0)

    def testAiteration(self):
        """Test the ability to iterate over the atoms"""
        atoms = [atom for atom in self.head[0]]
        self.assertEqual(len(atoms), self.Natoms)

    def Atomaccesstest(self):
        # Should raise AttributeError
        return self.atom.nosuchname

    def testAattributes(self):
        """Get the values of some properties"""
        self.assertRaises(AttributeError, self.Atomaccesstest)
        self.assertAlmostEqual(self.atom.coords[0], -0.0691, 4)

    def testAstringrepr(self):
        """Test the string representation of the Atom"""
        test = "Atom: 8 (-0.07 5.24 0.03)"
        self.assertEqual(str(self.atom), test)

    def testSMARTS(self):
        """Searching for ethyl groups in triethylamine"""
        mol = self.toolkit.readstring("smi", "CCN(CC)CC")
        smarts = self.toolkit.Smarts("[#6][#6]")
        ans = smarts.findall(mol)
        self.assertEqual(len(ans), 3)

class TestPybel(TestToolkit):
    toolkit = pybel
    tanimotoresult = 1/3.
    Ndescs = 3
    Natoms = 15
    tpsaname = "TPSA"
    Nbits = 3

    def testFP_FP3(self):
        "Checking the results from FP3"
        fps = [x.calcfp("FP3") for x in self.mols]
        self.assertEqual(fps[0] | fps[1], 0.)

    def testunitcell(self):
        """Testing unit cell access"""
        mol = self.toolkit.readfile("cif", "hashizume.cif").next()
        cell = mol.unitcell
        self.assertAlmostEqual(cell.GetAlpha(), 92.9, 1)

    def testMDcomment(self):
        """Mess about with the comment field"""
        data = self.head[0].data
        self.assertEqual('Comment' in data, True)
        self.assertEqual(data['Comment'], 'CORINA 2.61 0041  25.10.2001')
        data['Comment'] = 'New comment'
        self.assertEqual(data['Comment'], 'New comment')

    def testRSconversiontoMOL2(self):
        """Convert to mol2"""
        as_mol2 = self.mols[0].write("mol2")
        test = """@<TRIPOS>MOLECULE
*****
 4 3 0 0 0
SMALL
GASTEIGER
Energy = 0

@<TRIPOS>ATOM
      1 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      2 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      3 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
      4 C           0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
     3     3     4    1
"""
        self.assertEqual(as_mol2, test)

    def testRSgetprops(self):
        """Get the values of the properties."""
        self.assertAlmostEqual(self.mols[0].exactmass, 58.078, 3)
        self.assertAlmostEqual(self.mols[0].molwt, 58.122, 3)
        self.assertEqual(len(self.mols[0].atoms), 4)
        self.assertRaises(AttributeError, self.RSaccesstest)

class TestRDKit(TestToolkit):
    toolkit = rdkit
    tanimotoresult = 1/3.
    Ndescs = 176
    Natoms = 9
    tpsaname = "TPSA"
    Nbits = 12

    def testRSconversiontoMOL(self):
        """No conversion to MOL file done"""
        pass

  
class TestCDK(TestToolkit):
    toolkit = cdk
    tanimotoresult = 0.375
    Ndescs = 143
    Natoms = 15
    tpsaname = "tpsa"
    Nbits = 4

    def testSMARTS(self):
        """No SMARTS testing done"""
        pass

if __name__=="__main__":
    testcases = [TestPybel, TestCDK, TestRDKit]
    # testcases = [TestCDK]
    # testcases = [TestPybel]
    #testcases = [TestRDKit]
    for testcase in testcases:
        print "\n\n\nTESTING %s\n%s\n\n" % (testcase.__name__, "== "*10)
        myunittest = unittest.defaultTestLoader.loadTestsFromTestCase(testcase)
        unittest.TextTestRunner(verbosity=2).run(myunittest)
