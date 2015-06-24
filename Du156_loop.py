#!/usr/bin/env python

from modeller import *
from modeller.automodel import *
from modeller.parallel import *
from modeller.scripts import complete_pdb

# Loop refinement of an existing model

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = '../Du156_2B4C_4NCO_4TVP_mod/'
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        return selection(self.residue_range('349:A', '375:A'))

m = MyLoop(env,
           inimodel='Du156.B99990005.pdb', # initial model of the target
           sequence='Du156')          # code of the target

m.loop.starting_model= 1           # index of the first loop model 
m.loop.ending_model  = 10          # index of the last loop model
m.loop.md_level = refine.very_fast # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models

m.make()

for i in range(1, 11):
    # read model file
    code = "Du156.BL%04d0001.pdb" % i
    file = "Du156%04d.profile" %i
    mdl = complete_pdb(env, code)
    s = selection(mdl)
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=file,
                  normalize_profile=True, smoothing_window=15)
