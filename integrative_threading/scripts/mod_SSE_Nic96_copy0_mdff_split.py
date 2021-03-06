################################
# Full spoke, 6  helices
# Nic96 threading
################################

import IMP
import IMP.core
import IMP.atom
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.threading
import IMP.algebra
import IMP.pmi
import IMP.pmi.tools
import random
import numpy as np
import math
import itertools
import pandas as pd
import sys
from collections import OrderedDict
import timeit


from SSEThread import *
from IMP.pmi.tools import get_restraint_set

def groupbyUnsorted(input, key=lambda x:x):
  yielded = set()
  keys = [ key(element) for element in input ]
  for i, wantedKey in enumerate(keys):
    if wantedKey not in yielded:
      yield (wantedKey,
          (input[j] for j in range(i, len(input)) if keys[j] == wantedKey))
    yielded.add(wantedKey)

def all_start_iterator(lengths_new, limits_new, last_residue, sign, residue_interval = 4):
    all_lists = []
    for i in range(len(lengths_new)):
        all_lists.append(np.arange(limits_new[i][0], limits_new[i][1], residue_interval))
    
    all_start = []
    for c in itertools.product(*all_lists):
        c_ext = list(c)+[last_residue+1]
        diff = [(c_ext[i+1]-c_ext[i]) for i in range(len(c_ext)-1)]
        eval_diff = [d >=(l+2) for d,l in zip(diff, lengths_new)]
        if all(eval_diff):
            all_start.append(c)
    
    r_all_start = np.array(sum(all_start, ())).reshape(-1,len(lengths_new))
    return r_all_start

# Enumerate all threads
def get_limits(elements, first_residue, last_residue):
    # elements: (chain, first_residue, lenght)
    
    limits = []
    ri = first_residue
    rf = last_residue
    lengths = [int(e[2]) for e in elements]
    for i, l in enumerate(lengths[0:-1]):
        limits.append([ri, rf - sum(lengths[(i+1):])-l])
        ri = sum(lengths[0:(i+1)])+1
    limits.append([ri, rf-lengths[-1]+2])

    #print('Lengths', lengths)

    return limits

def eval_thread(SEs_ordered, SEMs_ordered, start_residues):
    for k, se_sel in enumerate(SEs_ordered):
        se_mover_sel = SEMs_ordered[k]
        se_mover_sel.zero_coordinates()
    
    for k, se_sel in enumerate(SEs_ordered):
        se_mover_sel = SEMs_ordered[k]
        l = se_sel.get_length()
        se_sel.set_start_res_key(start_residues[k])
        se_sel.set_length_key(int(l))
        se_mover_sel.transform_coordinates()

    cr_eval = cr.evaluate()
    print(cr_eval)
    if cr_eval < 5000:
        return cr_eval, xl.evaluate(), ppr_Nic96.evaluate(), ppr_Nup53.evaluate(), ppr_Nup59.evaluate()
    else:
        return 50001, np.nan, np.nan, np.nan, np.nan

def eval_threading_one(SEM, SEs_ordered,start_res):
    scores = []
    for st in start_res:
        for k, se_sel in enumerate(SEs_ordered):
            se_mover_sel = SEM[k]
        cre = 0
        cre = cr.evaluate()
        if cre > 5000.0:
            continue
        else:
            sc = list(st) + [cre] + [0]*3
            scores.append(sc)
            
    return scores

def eval_threading_sym(Tup, SE_ordered, start_residues, file_out):
    # Evaluate threading with 2-fold symmetry (2 copies of Nup53 and Nup59 and 4 copies Nic96)
    scores = []
    k = 0
    assign = [(a,b) for a,b,c in Tup]
    perm = [c for a,b,c in Tup]
    prots = [a for a,b,c in Tup]
    copies = [b for a,b,c in Tup]
    se_movers_sel = [SE_movers[(f'seq_chain_{prot}',copy)]['main'][k] for (prot,copy,k) in zip(prots,copies,perm)]
    se_movers_sym = [SE_movers[(f'seq_chain_{prot}',copy)]['sym'][k] for (prot,copy,k) in zip(prots,copies,perm)]
    print('Number of threading solutions',len(start_residues))
    for i, st in enumerate(start_residues):
            for k, se_sel in enumerate(SE_ordered):
                if i==0:
                    prot = Tup[k][0]
                    copy = Tup[k][1]
                    if prot == 'Nic96':
                        copy_sym = copy + 2
                    else:
                        copy_sym = copy + 1
                    l = int(se_sel[0].get_length())
                    se_sel[0].set_start_res_key(st[k])
                    se_sel[0].set_length_key(l)
                    se_sel[0].set_chain_key(f'seq_chain_{prot}_{copy}')
                    se_sel[1].set_start_res_key(st[k])
                    se_sel[1].set_length_key(l)
                    se_sel[1].set_chain_key(f'seq_chain_{prot}_{copy_sym}')
                    se_movers_sel[k].transform_coordinates()
                    se_movers_sym[k].transform_coordinates()
                    
                else:
                    if st[k] != start_residues[i-1][k]:
                        prot = Tup[k][0]
                        copy = Tup[k][1]
                        if prot == 'Nic96':
                            copy_sym = copy + 2
                        else:
                            copy_sym = copy + 1
                        se_movers_sel[k].zero_coordinates()
                        se_movers_sym[k].zero_coordinates()
                        prot = Tup[k][0]
                        copy = Tup[k][1] 
                        l = int(se_sel[0].get_length())
                        se_sel[0].set_start_res_key(st[k])
                        se_sel[0].set_length_key(l)
                        se_sel[0].set_chain_key(f'seq_chain_{prot}_{copy}')
                        se_sel[1].set_start_res_key(st[k])
                        se_sel[1].set_length_key(l)
                        se_sel[1].set_chain_key(f'seq_chain_{prot}_{copy_sym}')
                        se_movers_sel[k].transform_coordinates()
                        se_movers_sym[k].transform_coordinates()
                        
            cr_eval = cr.evaluate()
            if cr_eval < 4000:
                aa = [r for r in st]+[cr_eval, xl.evaluate(), ppr_Nic96.evaluate()]
                scores.append(aa)
            else:
                aa = list(st)+[cr_eval, np.nan, np.nan]
                scores.append(aa)

            if i%1000 == 0:
                print(i)
                write_scores_wcopies(assign, perm, scores, file_out)
                scores = []
    write_scores_wcopies(assign, perm, scores, file_out)
    return scores

def eval_threading(Tup, SE_ordered, start_residues, file_out):
    # Evaluate threading
    zero_coordinates_SEM(SE_movers[('seq_chain_Nic96',0)]['main'])  
  
    #all_polarities = list(itertools.product([-1, 1],repeat=len(SE_ordered)))
    all_polarities = [[1,1,1,1,1,1,1]]
    print('polarities',all_polarities, len(all_polarities), len(SE_ordered))
   
    scores = []
    k = 0
    assign = [(a,b) for a,b,c in Tup]
    perm = [c for a,b,c in Tup]
    prots = [a for a,b,c in Tup]
    copies = [b for a,b,c in Tup]
    se_movers_sel = [SE_movers[(f'seq_chain_{prot}',copy)]['main'][k] for (prot,copy,k) in zip(prots,copies,perm)]
    print('Number of threading solutions',len(start_residues))
    print('Tup', Tup)
    for i, st in enumerate(start_residues):
      for j, pol in enumerate(all_polarities):
        for k, se_sel in enumerate(SE_ordered):
            prot = Tup[k][0]
            copy = Tup[k][1]
            se_movers_sel[k].zero_coordinates()
            l = int(se_sel.get_length())
            se_sel.set_start_res_key(st[k])
            se_sel.set_length_key(l)
            if se_sel.get_polarity() != pol[k]:
              #print(j, k, pol[k])
              se_sel.flip_polarity_key()
            se_sel.set_chain_key(f'seq_chain_{prot}_{copy}')
            se_movers_sel[k].transform_coordinates()
            
        # Evaluate CR first
        cr_eval = cr.evaluate()
        if cr_eval < 80:
          print(st, cr_eval) 
          aa = [r for r in st]+[p for p in pol]+[cr_eval, xl.evaluate(), ppr_Nic96.evaluate()]
          scores.append(aa)
        else:
          aa = [r for r in st]+[p for p in pol]+[cr_eval, np.nan, np.nan]
          scores.append(aa)
        if i%100 == 0:
          write_scores_wcopies(assign, perm, scores, file_out)
          scores = []
        
    return scores
            
def zero_coordinates_SEM(SEM_ordered, SEM_sym = []):
    for sem in SEM_ordered+SEM_sym:
        sem.zero_coordinates()
    
def write_scores_wcopies(assign, perm, scores, file_out):
    assign_prot = (a[0] for a in assign)
    assign_copy = [a[1] for a in assign]
    p_prot = np.array(list(assign_prot) * len(scores)).reshape(-1, len(assign))
    p_perm = np.array(list(perm) * len(scores)).reshape(-1, len(perm))
    scores_concat = np.hstack((p_prot, p_perm, scores))
    np.savetxt(file_out, scores_concat, fmt='%s')

#########################
# Input information
#########################
print(len(sys.argv))
if len(sys.argv)<5:
  print('USAGE: python mod_SSE_Nic96_one_long_copy1.py permutation_selection[1-16] permutation_SE long_helix_pick[0-3]')
else:
  opt = int(sys.argv[1])
  sel_2 = int(sys.argv[2])
  resi_5 = int(sys.argv[3])
  resi_6 = int(sys.argv[4])
  print('permutation_selection[1-16] permutation_SE long_helix_pick[0-3]', opt, sel_2)

include_connectivity_restraint = True
include_distance_restraint = True
include_psipred_restraint  = True

############################
# Representation
############################
mdl = IMP.Model()

components = ['Nup157','Nup170','Nup53','Nup59','Nsp1','Nup57','Nup49','Nup188','Nup192' ]
colors = ['pink','orchid','gold','yellow','light green','green','sea green','plum','purple']
chains = [['Z','1'],['Y','0'],['U','W'],['V','X'],['A','D','G','J'],['B','E','H','K'],['C','F','I','L'],['N','P'],['M','O']]

beadsize=1

data_dir = '../../data/'
pdb_IR = f'{data_dir}/NPC_IR_fin-mdff_allcopies_orie.pdb'

# Setup System and add a State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()
mols = []

#IR -- single pdb
for i, (prot, chains_sel) in enumerate(zip(components,chains)):
    print('------- PMI: setting up',components[i])
    for k, ch in enumerate(chains_sel):
        if k == 0:
            seqs = IMP.pmi.topology.Sequences(f'{data_dir}/{components[i]}.fasta')
    
            mol = st.create_molecule(components[i],
                                     sequence=seqs[components[i]],
                                     chain_id=ch)
        
            atomic = mol.add_structure(pdb_IR,
                                       chain_id=ch,
                                       offset=0)
            mol.add_representation(atomic,                
                                   resolutions=[1],
                                   color=colors[i])
        else:
            print('Setting up clone:', prot, k)
            clon = mol.create_copy(chain_id = chains[i][k])
            atomic = clon.add_structure(pdb_IR,
                                        chain_id = chains[i][k])
            clon.add_representation(atomic,
                                    resolutions = [1],
                                    color=colors[i])


root_hier = s.build()

# Set all to optimized
for p in IMP.atom.Selection(root_hier).get_selected_particles():
    IMP.core.XYZ(p).set_coordinates_are_optimized(True)


# Write rmf to check that everthing matches
out = IMP.pmi.output.Output()
out.init_rmf("all_ini_fixed.rmf3", [root_hier])
out.write_rmf("all_ini_fixed.rmf3")
out.close_rmf("all_ini_fixed.rmf3")

############################
# Define structure elements
############################
seq_Nic96 = IMP.pmi.topology.Sequences(f'{data_dir}/Nic96.fasta')

seqs = {}
seqs['Nic96.0'] = seq_Nic96['Nic96']
seqs['Nic96.1'] = seq_Nic96['Nic96']

elements = [('Q',12,27-11),
            ('Q',36,40-35),
            ('Q',46,61-45),
            ('Q',65,77-64),
            ('Q',81,96-80),
            ('Q',134,157-133),
            ('Q',161,180-160)]

sec_struct = ['H','S','H','H','H','H','H']

SEs_chain = {'Q':[]}
for i, e in enumerate(elements):
  SEs_chain[e[0]].append(i)

print(SEs_chain)

print('lenghts', [e[2] for e in elements])

# Sequence to thread
chains = OrderedDict()
chains['Nic96'] = ['A',4,205]

SSE = SSEThread(root_hier = root_hier, sequences = seqs)
SEs_all = SSE.generate_structure_elements_from_hier(hier = root_hier, pdb=pdb_IR, elements= elements, sec_struct = sec_struct, stride = None)

print('all_SEs', SEs_all)

# Add fixed structure to sequence 

fixed_Nic96 = [('Q',206,839), ('R',206,839)]
SSE.add_coordinates_to_sequence_chain(hier = root_hier, prot = 'Nic96', pdb = pdb_IR, fixed=fixed_Nic96)

# Check hierarchy
#print(IMP.atom.show_with_representations(root_hier))

states = IMP.atom.get_by_type(root_hier, IMP.atom.STATE_TYPE)
mols = IMP.atom.get_by_type(root_hier, IMP.atom.MOLECULE_TYPE)
print('mols',mols)

seq_chains = SSE.get_sequence_hierarchy()
print('seq_chains:', seq_chains)

out.init_rmf("all_ini_Nic96.rmf3", [root_hier])
out.write_rmf("all_ini_Nic96.rmf3")
out.close_rmf("all_ini_Nic96.rmf3")

############################
# Setup for selected option
############################

# Three levels
# 1. Pick chains for SE - here always 0 
# 2. Pick chain for threading
# 3. Pick residues

chain_perms_sel = [['Q','Q','Q','Q','Q','Q','Q']]
print('Number of permutations:', len(chain_perms_sel),chain_perms_sel)

all_options = {}
prots = (('Nic96', 0))
for a, p in enumerate(chain_perms_sel):
    all_options[a] = (p,('Nic96', 0))
    print(a, p)
print(all_options, len(all_options))

print('selected option', all_options[opt])
# selected option, pick SE and generate labels
chains_label = ''.join(all_options[opt][0])
prot_copy = all_options[opt][1][0]

prots = (all_options[opt][1],)*7

SEs = {}
elements_single = []
for i, chain in enumerate(all_options[opt][0]):
  SE_index = SEs_chain[chain][i]
  SEs[i] = SEs_all[SE_index]
  elements_single.append(elements[SE_index])

print(elements_single)
print(chains_label, SEs)

# Now all permutations 
sel_assign = {}
perms_key = {}
perms_key['Q'] = [(0,1,2,3,4,5,6)]
    
sel_assign[tuple(all_options[opt][0])] = perms_key
print(perms_key, all_options)
    

sel_assign[tuple(all_options[opt][0])] = perms_key

# Now all products
for assign, perms in sel_assign.items():
    print('assign', assign)
    print(len(list(itertools.product(*perms.values()))),list(itertools.product(*perms.values())))

print(sel_assign[tuple(all_options[opt][0])].keys())

############################
# Setup Movers
############################
se_movers = []
seq_chains = SSE.get_sequence_hierarchy()

copy_ch = all_options[opt][1][1]
print('copy', copy_ch, all_options[opt])

chain = [ch for ch in seq_chains if IMP.atom.Copy(ch).get_copy_index()==copy_ch][0]

SE_movers = {}
SE_movers[('seq_chain_Nic96',0)] = {'main':[]}
for i, se in SEs.items():
  semi = SEMover(SSE, k, chid=copy_ch)
  print('index', i, se.get_particle_index())
  sem = IMP.threading.StructureElementMover(mdl, se.get_particle_index(), chain.get_particle())
  SE_movers[('seq_chain_Nic96',copy_ch)]['main'].append(sem)
  sem.zero_coordinates()

print(SE_movers)

############################
# Initialize
############################

prot = prots[0][0]
copy = prots[0][1]

st = [12,32,46,65,81,134,161]
perm = [0,1,2,3,4,5,6]
se_movers_sel = [SE_movers[(f'seq_chain_{prot}',copy)]['main'][k] for k in perm]


for k, se_sel in SEs.items():
  l = int(se_sel.get_length())
  se_sel.set_start_res_key(st[k])
  se_sel.set_length_key(l)
  se_sel.set_chain_key(f'seq_chain_{prot}_{copy}')
  se_movers_sel[k].transform_coordinates()

print(IMP.atom.show_with_representations(root_hier))

print('lenghts', [e[2] for e in elements])

############################
# Connectivy restraint
############################
restraints = []

Nic96_C_fixed = IMP.atom.Selection(root_hier, chain_id='SNic96', residue_index=206).get_selected_particles()
print('fixed', Nic96_C_fixed)
cr = ChainConnectivityRestraint(root_hier,
                                SEs.values(),
                                residues_C = Nic96_C_fixed)
cr.add_to_model()
restraints.append(cr)


############################
# Cross-linking restraint
############################
print('Setting up cross-linking restraint ...')
chain_mapping = {'Nic96':'S0'}
#dir = '/Users/iecheverria/Dropbox/UCSF/yeast_npc/structural_coverage'
xl = StructuralElementCrossLinksRestraint(root_hier,
                                          f'{data_dir}/XLs_all_2020_sel_Nic96_Nup53_Nup59_Nup100_Nup116_Nup145N.csv',
                                          chain_mapping)
xl.add_to_model()
restraints.append(xl)
xl_score = xl.evaluate()
print('xl score', xl.evaluate())


###############################
# Setup SS restraint
###############################
print('Setting up Secondary Structure restraint ...')
ppr_Nic96 = SecondaryStructureRestraint(hierarchy=root_hier,
                                        structural_elements=SEs,
                                        prot = 'Nic96', 
                                        ss_file = f'{data_dir}/RaptorX_Nic96.txt',
                                        input_format = 'RaptorX',
                                        label='SSP_Nic96')

ppr_Nic96.set_weight(1.0)
ppr_Nic96.add_to_model()
restraints.append(ppr_Nic96)
def run_ppr_Nic96():
    return ppr_Nic96.evaluate()

print('SS score Nic96:', ppr_Nic96.evaluate())

###########################
# Time everything
###########################
t_cr = timeit.timeit(run_cr,number=1)
t_xl = timeit.timeit(run_xl,number=1)
t_Nic96 = timeit.timeit(run_ppr_Nic96,number=1)

print('cr, xl, Nip96',t_cr,t_xl,t_Nic96)

###########################
# For testing
###########################
#12
out_perm = open('test.dat','w')
perm = (0,1,2,4,3,5,6)

SEs_ordered = [SEs[i] for i in perm]
print('SEs_ordered', SEs_ordered)

start_res = [[10., 31., 48.0, 66.0, 93.0, 119.0, 142.0],[10., 31.0, 48.0, 66.0, 93.0, 119.0, 142.0], [10., 31.0, 50.0, 66.0, 93.0, 119.0, 142.0], [10., 31.0, 50.0, 66.0, 93.0, 119.0, 148.0]]

Tup = [(a,b,c) for (a,b),c in zip(prots,perm)]
print('Tup', Tup)
scores_all = eval_threading(Tup, SEs_ordered, start_res, out_perm)
out_perm.close()

print(IMP.atom.show_with_representations(root_hier))



print(xl.evaluate())
print(all_options[opt], opt, np.array(scores_all))
zero_coordinates_SEM(SE_movers[('seq_chain_Nic96',0)]['main'])
###########################
# Scoring function
############################
RS = get_restraint_set(mdl)
sf = IMP.core.RestraintsScoringFunction(RS)

############################
# All possible protein
# assignments
# Filter assignments to have
# the same protein for SE 2
# and 3
############################

first_residue = chains['Nic96'][1]
last_residue = chains['Nic96'][2]

n_models = 0

print('sel_assign', sel_assign)
for assign, perms in sel_assign.items():
    # Simple dictionary to keep chain counts
    D = {'Q' : 0, 'R' : 0, 'S':0}
    nassign = []
    for a in assign:
      nassign.append(D[a])
      D[a]+=1
    zero_coordinates_SEM(SE_movers[('seq_chain_Nic96',0)]['main'])
    n_models_it = 0
    print(len(list(itertools.product(*perms.values()))))
    for ii, iters in enumerate(itertools.product(*perms.values())):
      if ii == sel_2:
        iters_d = {k:d for k,d in zip(perms.keys(),iters)}
        permutation = [iters_d[ch][i] for i, ch in zip(nassign, assign)]
        print('permutation', permutation)
        perm_label = ''.join([str(p) for p in permutation])
        print(assign, permutation, chains_label)
        out_assign = open(f'enumeration_perm_{chains_label}_copy{copy_ch}_{perm_label}_resi{resi_5}_{resi_6}.csv','w')
        print('iters', iters, permutation)
        lengths_new = [elements_single[i][2] for i in permutation]
        print('lenghts', lengths_new)
        elements_new = [elements_single[i] for i in permutation]
        limits_new = get_limits(elements_new, first_residue, last_residue)
        limits_new[0][0]=10
        limits_new[0][1]=20
        limits_new[1][0]=29
        limits_new[1][1]=42
        limits_new[2][0]=45
        limits_new[2][1]=58
        limits_new[3][0]=62
        limits_new[3][1]=80
        limits_new[4][0]=80
        limits_new[4][1]=96
        limits_new[5][0]=resi_5
        limits_new[5][1]=resi_5+1
        limits_new[6][0]=resi_6
        limits_new[6][1]=resi_6+1
        print(limits_new)
        sign = list(np.sign(np.diff(permutation)))+[1]
        start_res = all_start_iterator(lengths_new, limits_new, last_residue, sign, residue_interval = 1)
        print('first_row', start_res[0:10,:])
        print('last_row', start_res[-1,:])
        print('shape1', np.shape(start_res))
        print([start_res[:,i+1]-start_res[:,i] for i in range(6)]) 
        print(np.array([start_res[:,i+1]-start_res[:,i]-lengths_new[i] for i in range(6)]).T)
        Tup = [(p[0], p[1], ii) for p, ii in zip(prots, permutation)]
        print('here', Tup, np.shape(start_res), assign, len(assign))
        SEs_ordered = [SEs[i] for i in permutation]
        # Finally, evaluate scores
        scores_all = eval_threading(Tup, SEs_ordered, start_res, out_assign)
        out_assign.close()
        
print('Completed!!!!') 

exit()


