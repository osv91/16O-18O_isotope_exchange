from rdkit import Chem
from rdkit.Chem.inchi import *

class fragmentation_tree(inchi, fragment_formulas):
    def __init__(self, inchi, fragment_formulas):
        self.mol = Chem.MolFromInchi(self.inchi)