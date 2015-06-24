#!/usr/bin/env python

from modeller import *
from modeller.automodel import automodel
from modeller.automodel import assess
from modeller.parallel import *
import sys

def run_modeller(knowns = ('2B4C', '4NCO'), seq = 'CAP45',
                  num_mod = 1):
    
    seq_dir = '../sequence_files/'
    pdb_dir = '../pdb_files/'
    
    j = job()
    j.append(local_slave())
    j.append(local_slave())
    j.append(local_slave())
    j.append(local_slave())
    j.append(local_slave())

    # Modeller environment
    env = environ()
    
    # Output
    log.none()
    
    #input dir
    env.io.atom_files_directory = [pdb_dir]
        
    #input file
    ali_file = seq_dir + 'align' + seq
    if isinstance(knowns, list):
        for i in range(len(knowns)):
            ali_file += '_' + knowns[i]
        ali_file += '.ali'
    else:
        ali_file += '_' + str(knowns) + '.ali'
    
    mod = automodel(env,
                    alnfile = ali_file,
                    knowns = knowns,
                    sequence = seq,
                    assess_methods = (assess.DOPE))

    mod.starting_model = 1
    mod.ending_model = num_mod

    mod.use_parallel_job(j)
    mod.make() 
    
    # Get list of all built models
    ok_models = [x for x in mod.outputs if x['failure'] is None]

    # Rank the models by DOPE score
    ok_models.sort(key = lambda mod: mod['DOPE score'])

    # Get top model
    m = ok_models[0]
    print("The best model {} with a dope of {}".format(m['name'], m['DOPE score']))

seq = sys.argv[1]
num_mod = int(sys.argv[2])

knowns = []
for i in range(3, len(sys.argv)):
    knowns.append(sys.argv[i])
    print(knowns)


run_modeller(knowns = knowns, seq = seq, num_mod = num_mod)
