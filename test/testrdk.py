import os
import unittest
import rdk as pybel

# For compatability with Python2.3
try:
    from sets import Set as set
except ImportError:
    pass

class Test_fingerprint(unittest.TestCase):
    """Test the Fingerprint class"""
    def setUp(self):
        self.mols = [pybel.readstring("smi", "CCCC"),
                     pybel.readstring("smi", "CCCN")]

    def testTanimoto(self):
        """Test the calculation of the Tanimoto coefficient"""
        fps = [x.calcfp() for x in self.mols]
        self.assertEqual(fps[0] | fps[1], 1/3.)
        fps = [x.calcfp("FP3") for x in self.mols]
        self.assertEqual(fps[0] | fps[1], 0.)
    
    def teststringrepr(self):
        """Test the string representation and corner cases."""
        self.assertRaises(ValueError, self.mols[0].calcfp, "Nosuchname")
        self.assertRaises(AttributeError, self.accesstest)
        self.assertEqual(str(self.mols[0].calcfp()),
"33554432, 0, 0, 0, 0, 0, 0, 2048, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,"
" 0, 0, 0, 0, 1073741824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2097152, 0, 131072,"
" 0, 0, 0, 2147483680, 0, 0, 0, 0, 0, 0, 0, 0, 64, 0, 0, 0, 8, 0, 2048, 0, 0, 0,"
" 0, 1, 0, 0")

    def accesstest(self):
        # Should raise AttributeError
        return self.mols[0].calcfp().nosuchname

    def testbits(self):
        """Test whether the bits are set correctly."""
        bits = [x.calcfp().bits for x in self.mols]
        self.assertEqual(bits[0],
                         [25, 235, 257, 734, 1141, 1201, 1317,
                          1343, 1606, 1731, 1803, 1952])
        bits = [set(x) for x in bits]
        # Calculate the Tanimoto coefficient the old-fashioned way
        tanimoto = len(bits[0] & bits[1]) / float(len(bits[0] | bits[1]))
        self.assertEqual(tanimoto, 1/3.)

class Test_readstring(unittest.TestCase):
    """Test the ability to read and write to a string"""
    def setUp(self):
        self.mol = pybel.readstring("smi", "CCCC")
        
    def accesstest(self):
        # Should raise AttributeError
        return self.mol.nosuchname

    def testformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, pybel.readstring, "noel", "jkjk")
    
    def testgetprops(self):
        """Get the values of the properties."""
        assert abs(self.mol.molwt - 58.121) < 0.003
        self.assertEqual(len(self.mol.atoms), 4)
        self.assertRaises(AttributeError, self.accesstest)

    def teststringrepr(self):
        """Test the string representation of a molecule"""
        self.assertEqual(str(self.mol).strip(), "CCCC")

class Test_readfile(unittest.TestCase):
    """Test the ability to read and write to a file"""
    def setUp(self):
        self.mols = [mol for mol in pybel.readfile("sdf", "head.sdf")]

    def testread(self):
        """Is the right number of molecules read from the file?"""
        self.assertEqual(len(self.mols), 2)

    def readerror(self):
        mol = pybel.readfile("sdf", "nosuchfile.sdf").next()

    def testmissingfile(self):
        """Test that reading from a non-existent file raises an error."""
        self.assertRaises(IOError, self.readerror)

    def formaterror(self):
        mol = pybel.readfile("noel", "head.sdf").next()
    
    def testformaterror(self):
        """Test that invalid formats raise an error"""
        self.assertRaises(ValueError, self.formaterror)

    def unitcellerror(self):
        unitcell = self.mols[0].unitcell
    
    def testunitcellerror(self):
        """Test that accessing the unitcell raises an error"""
        self.assertRaises(AttributeError, self.unitcellerror)

    def testconversion(self):
        """Convert to smiles"""
        as_smi = [mol.write("smi").split("\t")[0] for mol in self.mols]
        test = ['CC1=CC(=O)C=CC1=O', 'c1ccc2c(c1)nc(SSc1nc3ccccc3s1)s2']
        self.assertEqual(as_smi, test)

    def test_singletofile(self):
        """Test the molecule.write() method"""
        mol = self.mols[0]
        mol.write("smi", "testoutput.txt")
        test = ['CC1=CC(=O)C=CC1=O\n']
        filecontents = open("testoutput.txt", "r").readlines()
        self.assertEqual(filecontents, test)
        self.assertRaises(IOError, mol.write, "smi", "testoutput.txt")
        os.remove("testoutput.txt")
        self.assertRaises(ValueError, mol.write, "noel", "testoutput.txt")
    
    def test_multipletofile(self):
        """Test the Outputfile class"""
        self.assertRaises(ValueError, pybel.Outputfile, "noel", "testoutput.txt")
        outputfile = pybel.Outputfile("smi", "testoutput.txt")
        for mol in self.mols:
            outputfile.write(mol)
        outputfile.close()
        self.assertRaises(IOError, outputfile.write, mol)
        self.assertRaises(IOError, pybel.Outputfile, "smi", "testoutput.txt")
        filecontents = open("testoutput.txt", "r").readlines()
        os.remove("testoutput.txt")
        test = ['SMILES Name \n',
                'CC1=CC(=O)C=CC1=O NSC 1\n',
                'c1ccc2c(c1)nc(SSc1nc3ccccc3s1)s2 NSC 2\n']
        self.assertEqual(filecontents, test)

    def desctest(self):
        # Should raise ValueError
        self.mols[0].calcdesc("BadDescName")

    def testdesc(self):
        """Test the descriptors"""
        desc = self.mols[0].calcdesc()
        self.assertEqual(len(desc), 176)
        self.assertAlmostEqual(desc['MolLogP'], 0.64, 2)
        self.assertRaises(ValueError, self.desctest)

class Test_data(unittest.TestCase):
    def setUp(self):
        self.mol = pybel.readfile("sdf", "head.sdf").next()
        self.data = self.mol.data

    def accesstest(self):
        # Should raise KeyError
        return self.data['noel']

    def testaccess(self):
        """Change the value of a field"""
        self.assertRaises(KeyError, self.accesstest)
        self.data['noel'] = 'testvalue'
        self.assertEqual(self.data['noel'], 'testvalue')
        newvalues = {'hey':'there', 'yo':1}
        self.data.update(newvalues)
        self.assertEqual(self.data['yo'], '1')
        self.assert_('there' in self.data.values())

    def testglobalaccess(self):
        """Check out the keys"""
        self.assert_(self.data.has_key('NSC'))
        self.assert_(not self.data.has_key('Noel'))
        self.assertEqual(len(self.data), 1)
        for key in self.data:
            self.assertEqual(key in ['NSC'], True)
        self.assertEqual(repr(self.data), "{'NSC': '1'}")

    def testdelete(self):
        """Delete some keys"""
        self.assert_(self.data.has_key('NSC'))
        del self.data['NSC']
        self.assert_(not self.data.has_key('NSC'))
        self.data.clear()
        self.assertEqual(len(self.data), 0)

class Test_atoms(unittest.TestCase):
    """Testing some of the atom code"""
    def setUp(self):
        self.mol = pybel.readfile("sdf", "head.sdf").next()
        self.atom = self.mol.atoms[0]

    def testiteration(self):
        """Test the ability to iterate over the atoms"""
        atoms = [atom for atom in self.mol]
        self.assertEqual(len(atoms), 9)
        self.mol.addh()
        atoms = [atom for atom in self.mol]
        self.assertEqual(len(atoms), 15)        

    def accesstest(self):
        # Should raise AttributeError
        return self.atom.nosuchname

    def testattributes(self):
        """Get the values of some properties"""
        self.assertRaises(AttributeError, self.accesstest)
        self.assert_(abs(self.atom.coords[0]-0.0021) < 0.0001)

    def teststringrepr(self):
        """Test the string representation of the Atom"""
        test = "Atom: 8 (0.0021, -0.0041, 0.0020)"
        self.assertEqual(str(self.atom), test)

class Test_smarts(unittest.TestCase):
    """Test the Smarts object"""
    def setUp(self):
        self.mol = pybel.readstring("smi", "CCN(CC)CC")

    def testmatching(self):
        """Searching for ethyl groups in triethylamine"""
        smarts = pybel.Smarts("[#6][#6]")
        ans = smarts.findall(self.mol)
        self.assertEqual(len(ans), 3)
    
if __name__=="__main__":
    unittest.main()
