from modeller import *
from modeller.automodel import *
import sys

# Override the 'special_restraints' and 'user_after_single_model' methods:

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        
    def special_patches(self, aln):
        self.rename_segments(segment_ids=('A'), renumber_residues=[260])

env = environ()

env.io.hetatm = True

a = MyModel(env, alnfile='aln_ig1_hhpred.pir', knowns=('5tvz'),sequence='ypom152_ig1',assess_methods=(assess.DOPE, assess.normalized_dope))

a.starting_model = 1
a.ending_model = 50
if '--test' in sys.argv: a.ending_model = 1

a.make()                           # do comparative modeling

print()
