class Molecule(object):
    """Superclass of cinfony molecules.

    Should not be instantiated directly.
    """

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.
        
        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        for atom in self.atoms:
            yield atom

    def __str__(self):
        return self.write("smi")

if __name__=="__main__": #pragma: no cover
    pass
    
