######################################################
# Full spoke, 4 long helices
# Nup53, Nup59, Nup100, Nup116, Nup145N  threading
######################################################

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
        eval_diff = [d >=(l-1) for d,l in zip(diff, lengths_new)]
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
        limits.append((ri, rf - sum(lengths[(i+1):])-l))
        ri = sum(lengths[0:(i+1)])+1
    limits.append((ri, rf-lengths[-1]+2))

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
        #se_sel.set_chain_key(chain_id[prots[k]])
        se_mover_sel.transform_coordinates()

    cr_eval = cr.evaluate()
    print(cr_eval)
    if cr_eval < 5000:
        return cr_eval, xl.evaluate(),  ppr_Nup53.evaluate(), ppr_Nup59.evaluate()
    else:
        return 50001, np.nan, np.nan, np.nan, np.nan

def eval_threading_one(SEM, SEs_ordered,start_res):
    scores = []
    for st in start_res:
        for k, se_sel in enumerate(SEs_ordered):
            se_mover_sel = SEM[k]
            #se_mover_sel.zero_coordinates()
            #l = se_sel.get_length()
            #se_sel.set_start_res_key(start_residues[k])
            #se_sel.set_length_key(int(l))
            #se_mover_sel.transform_coordinates()
        cre = 0
        cre = cr.evaluate()
        if cre > 5000.0:
            continue
        else:
            sc = list(st) + [cre] + [0]*3
            #sc = start_res + [cre] + [r.evaluate() for r in other_restraints]
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
    se_movers_sel = [SE_movers[(f'seq_chain_{prot}',copy)][k] for (prot,copy,k) in zip(prots,copies,perm)]
    #se_movers_sym = [SE_movers[(f'seq_chain_{prot}',copy)]['sym'][k] for (prot,copy,k) in zip(prots,copies,perm)]
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
                aa = [r for r in st]+[cr_eval, np.nan, np.nan]
                scores.append(aa)

            if i%1000 == 0:
                print('i', i)
                write_scores_wcopies(assign, perm, scores, file_out)
                scores = []
    write_scores_wcopies(assign, perm, scores, file_out)
    return scores

def eval_threading(Tup, SE_ordered, start_residues, file_out):
    # Evaluate threading

    all_polarities = [p for p in itertools.product([-1, 1],repeat=2)]

    scores = []
    k = 0
    assign = [(a,b) for a,b,c in Tup]
    perm = [c for a,b,c in Tup]
    prots = [a for a,b,c in Tup]
    copies = [b for a,b,c in Tup]
    print('Tup, polarities', Tup, all_polarities)
    se_movers_sel = [SE_movers[(f'seq_chain_{prot}',copy)][k] for (prot,copy,k) in zip(prots,copies,perm)]
    print('Number of threading solutions',len(start_residues))
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
              se_sel.flip_polarity_key()
            se_sel.set_chain_key(f'seq_chain_{prot}_{copy}')
            se_movers_sel[k].transform_coordinates()
        # Evaluate CR first
        cr_eval = cr.evaluate()
        if cr_eval < 80:
          aa = [r for r in st]+[p for p in pol]+[cr_eval, xl.evaluate()]+[ppr.evaluate() for ppr in ppr_restraints] 
          scores.append(aa)
        else:
          aa = [r for r in st]+[p for p in pol]+[cr_eval, np.nan]+ [np.nan for i in range(len(ppr_restraints))]

          scores.append(aa)
      if i%100 == 0:
        write_scores_wcopies(assign, perm, scores, file_out)
        scores = []
    # Zero coordinates
    for prot, copy, pp in  Tup:
      zero_coordinates_SEM(SE_movers[(f'seq_chain_{prot}',copy)])
    write_scores_wcopies(assign, perm, scores, file_out)
    return scores

def zero_coordinates_SEM(SEM_ordered, SEM_sym = []):
    for sem in SEM_ordered+SEM_sym:
        sem.zero_coordinates()
    
def write_scores_wcopies(assign, perm, scores, file_out):
    assign_prot = (a[0] for a in assign)
    assign_copy = [a[1] for a in assign]
    p_prot = np.array(list(assign_prot) * len(scores)).reshape(-1, len(assign))
    p_copy = np.array(list(assign_copy) * len(scores)).reshape(-1, len(assign))
    p_perm = np.array(list(perm) * len(scores)).reshape(-1, len(perm))
    scores_concat = np.hstack((p_prot,p_copy, p_perm, scores))
    np.savetxt(file_out, scores_concat, fmt='%s')

def clean_XLs_database(file_in, prots_in, sel_int):
    file_name = file_in.split('.')[0]
    out = open(f'{file_name}_updated_{sel_int}.csv','w')
    xls_eff = 0
    xls_ign = 0
    for line in open(file_in):
        vals = line.split(',')
        if vals[0] == 'Protein1':
            out.write(line)
        else:
            prot1 = vals[0]
            prot2 = vals[2]
            r1 = int(vals[1])
            r2 = int(vals[3])

            ps1 = IMP.atom.Selection(root_hier,
                                     molecule=prot1,
                                     residue_index=r1).get_selected_particles()
            ps2 = IMP.atom.Selection(root_hier,
                                     molecule=prot2,
                                     residue_index=r2).get_selected_particles()

            if len(ps1) == 0 and prot1 in prots_in:
              out.write(line)
              xls_eff += 1
              continue
            elif len(ps2) == 0 and prot2 in prots_in:
              out.write(line)
              xls_eff += 1
    out.close()
    print('XLs effective', xls_eff )

#########################
# Input information
#########################
pp = sys.argv[1]
sel_2 = int(sys.argv[2])

include_connectivity_restraint = True
include_distance_restraint = True
include_psipred_restraint  = True

############################
# Representation
############################
mdl = IMP.Model()

components = ['Nup157','Nup170','Nsp1','Nup57','Nup49','Nup188','Nup192','Nic96']
colors = ['pink','orchid','light green','green','sea green','plum','purple','light blue']
chains = [['Z','1'],['Y','0'],['A','D','G','J'],['B','E','H','K'],['C','F','I','L'],['N','P'],['M','O'],['Q','R','S','T']]

beadsize=1

data_dir = 'data/'
data_dir = '/wynton/home/sali/ignacia/NPC/SSE/SSE_Nup192_helices_mdff/data'
pdb_IR = f'{data_dir}/NPC_IR_fin-mdff_allcopies_orie_renumbered.pdb'
pdb_cytoplasm = f'{data_dir}/Cytoplasm_rebuilt.pdb'
pdb_nucleoplasm = f'{data_dir}/Nucleoplasm_rebuilt.pdb'
pdb_orphan = f'{data_dir}/main_ISR_OSD_Nup192_orphan_helix.pdb'


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
seq_Nup53 = IMP.pmi.topology.Sequences(f'{data_dir}/Nup53.fasta')
seq_Nup59 = IMP.pmi.topology.Sequences(f'{data_dir}/Nup59.fasta')
seq_Nup100 = IMP.pmi.topology.Sequences(f'{data_dir}/Nup100.fasta')
seq_Nup116 = IMP.pmi.topology.Sequences(f'{data_dir}/Nup116.fasta')
seq_Nup145N = IMP.pmi.topology.Sequences(f'{data_dir}/Nup145N.fasta')

seqs = {}
seqs['Nup53.0'] = seq_Nup53['Nup53']
seqs['Nup53.1'] = seq_Nup53['Nup53']
seqs['Nup59.0'] = seq_Nup59['Nup59']
seqs['Nup59.1'] = seq_Nup59['Nup59']

seqs['Nup100.0'] = seq_Nup100['Nup100']
seqs['Nup100.1'] = seq_Nup100['Nup100']

seqs['Nup116.0'] = seq_Nup116['Nup116']
seqs['Nup116.1'] = seq_Nup116['Nup116']

seqs['Nup145N.0'] = seq_Nup145N['Nup145N']
seqs['Nup145N.1'] = seq_Nup145N['Nup145N']


elements = [('T',1,36),('R',1,36)]
sec_struct = ['H','H']

# Sequence to thread
chains = {}
chains['Nup53'] = ['B', 2,247-36]
chains['Nup59'] = ['B', 2,265-36]
chains['Nup100'] = ['A', 551,815-36]
chains['Nup116'] = ['A', 751,965-36]
chains['Nup145N'] = ['A', 201,458-36]

SSE = SSEThread(root_hier = root_hier, sequences = seqs)
SEs = SSE.generate_structure_elements_from_hier(hier = root_hier,
                                                pdb=pdb_orphan,
                                                elements= elements,
                                                sec_struct = sec_struct,
                                                stride = None)


# Add fixed structure to sequence
fixed_Nup53 = [('U',248,360),('W',248,360)]
SSE.add_coordinates_to_sequence_chain(hier = root_hier, prot = 'Nup53', pdb = pdb_IR, fixed=fixed_Nup53)

fixed_Nup59 = [('V',266,402),('X',266,402)]
SSE.add_coordinates_to_sequence_chain(hier = root_hier, prot = 'Nup59', pdb = pdb_IR, fixed=fixed_Nup59)

fixed_Nup100 = [('C',816,959),('D',816,959)]
SSE.add_coordinates_to_sequence_chain(hier = root_hier, prot = 'Nup100', pdb = pdb_cytoplasm, fixed=fixed_Nup100)

fixed_Nup116 = [('A',966,1111),('B',966,1111)]
SSE.add_coordinates_to_sequence_chain(hier = root_hier, prot = 'Nup116', pdb = pdb_cytoplasm, fixed=fixed_Nup116)

fixed_Nup145N = [('A',459,605),('B',459,605)]
SSE.add_coordinates_to_sequence_chain(hier = root_hier, prot = 'Nup145N', pdb = pdb_nucleoplasm, fixed=fixed_Nup145N)

# Check hierarchy
#print(IMP.atom.show_with_representations(root_hier))

states = IMP.atom.get_by_type(root_hier, IMP.atom.STATE_TYPE)
mols = IMP.atom.get_by_type(root_hier, IMP.atom.MOLECULE_TYPE)
print('mols',mols)

seq_chains = SSE.get_sequence_hierarchy()
print(seq_chains)

out.init_rmf("all_ini.rmf3", [root_hier])
out.write_rmf("all_ini.rmf3")
out.close_rmf("all_ini.rmf3")

############################
# Setup Movers
############################
se_movers = []
seq_chains = SSE.get_sequence_hierarchy()
print(seq_chains)


# Need to generate one SEMover for each chain/SE 
SE_movers = {k:[] for k in chains.keys()}

print('chains', seq_chains)
print('SEs', SEs)

SE_movers = {}
SE_movers[('seq_chain_Nup53',0)] = []
SE_movers[('seq_chain_Nup53',1)] = []
SE_movers[('seq_chain_Nup59',0)] = []
SE_movers[('seq_chain_Nup59',1)] = []
SE_movers[('seq_chain_Nup100',0)] = []
SE_movers[('seq_chain_Nup100',1)] = []
SE_movers[('seq_chain_Nup116',0)] = []
SE_movers[('seq_chain_Nup116',1)] = []
SE_movers[('seq_chain_Nup145N',0)] = []
SE_movers[('seq_chain_Nup145N',1)] = []

for i, ch in enumerate(seq_chains):
  copy_ch =  IMP.atom.Copy(ch).get_copy_index()
  for k, se in SEs.items():
    semi = SEMover(SSE, k, chid=copy_ch)
    sem = IMP.threading.StructureElementMover(mdl, se.get_particle_index(), ch.get_particle())
    SE_movers[(ch.get_name(),copy_ch)].append(sem)         
    sem.zero_coordinates()
                
print('SE_movers', SE_movers)
############################
# All possible protein
# assignments
# Filter assignments to have
# the same protein for SE 2
# and 3
############################

all_assign = {}

prots_sel = [('Nup53',0),('Nup53',1),
             ('Nup59',0),('Nup59',1),
             ('Nup100',0),('Nup100',1),
             ('Nup116',0),('Nup116',1),
             ('Nup145N',0),('Nup145N',1)]

k = 0
for i,prod in enumerate(itertools.product(prots_sel,repeat = 2)):
  perms_key = {}
  if len(set(prod))==1:
    continue
  else:
    # Eliminate the ones that have only one Nup53, Nup59
    n_Nup53 = np.sum([1 for p in prod if p[0]=='Nup53'])
    n_Nup59 = np.sum([1 for p in prod if p[0]=='Nup59'])
    if n_Nup53 == 1 or n_Nup59==1:
      continue
    else:
      print(prod)
      k += 1
      # get positions per protein  
      for key, subiter in groupbyUnsorted(list(enumerate(prod)),lambda x: x[1]):
        info = list(subiter)
        pos_key = [x[0] for x in info]
        perms_key[key] = [p for p in itertools.permutations(pos_key, len(pos_key))] 
    all_assign[prod] = perms_key

print('Number of assignments:', len(all_assign))

tot_n_it = 0
for k,v in all_assign.items():
  
  tot_n_it += len(v)
print(len(all_assign), tot_n_it)

sel_assign = {}
print('all_assign', len(all_assign))
for i, (k, v) in enumerate(all_assign.items()):
    all_p = list(itertools.product(*v.values()))
    print('----', i, k, len(v), (len(all_p)))
    if i == int(pp):
        sel_assign[k] = v
        if len(all_p) < (sel_2+1):
          exit()
print('sel_assign', sel_assign)
exit()
############################
# Connectivy restraint
############################
restraints = []

Nup53_C_fixed = IMP.atom.Selection(root_hier,chain_id='SNup53', residue_index=248).get_selected_particles()
Nup59_C_fixed = IMP.atom.Selection(root_hier,chain_id='SNup59', residue_index=266).get_selected_particles()
Nup100_C_fixed = IMP.atom.Selection(root_hier,chain_id='SNup100', residue_index=816).get_selected_particles()
Nup116_C_fixed = IMP.atom.Selection(root_hier,chain_id='SNup116', residue_index=966).get_selected_particles()
Nup145N_C_fixed = IMP.atom.Selection(root_hier,chain_id='SNup145N', residue_index=459).get_selected_particles()

residues_C = Nup53_C_fixed + Nup59_C_fixed + Nup100_C_fixed + Nup116_C_fixed + Nup145N_C_fixed

print(residues_C, len(residues_C))

cr = ChainConnectivityRestraint(root_hier,
                                SEs.values(),
                                residues_C = residues_C)
cr.add_to_model()
restraints.append(cr)

def run_cr():
    return cr.evaluate()
print('--------')
#print('Connectivity score:', cr.evaluate())
#exit()
############################
# Cross-linking restraint
############################
sel_proteins = list(set([p[0] for k in sel_assign.keys() for p in k]))
print('Setting up cross-linking restraint ...')
clean_XLs_database(data_dir+"/XLs_all_2020.csv",
                   sel_proteins, pp)

print(sel_proteins)

xl = StructuralElementCrossLinksRestraint(root_hier,
                                          f'{data_dir}/XLs_all_2020_updated_{pp}.csv')

xl.add_to_model()
restraints.append(xl)
xl_score = xl.evaluate()
print('xl score', xl.evaluate())

def run_xl():
    return xl.evaluate()

###############################
# Setup SS restraint
###############################
print('Setting up Secondary Structure restraint ...')

print(sel_assign)

all_proteins = ['Nup53','Nup59','Nup100','Nup116','Nup145N']
sel_proteins = list(set([p[0] for k in sel_assign.keys() for p in k]))
ppr_restraints = []

for prot in all_proteins:
  if prot in sel_proteins:
    ppr = SecondaryStructureRestraint(hierarchy=root_hier,
                                      structural_elements=SEs,
                                      prot = prot, 
                                      ss_file = f'{data_dir}/RaptorX_{prot}.txt',
                                      input_format = 'RaptorX',
                                      label=f'SSP_{prot}')

    ppr.set_weight(1.0)
    ppr_restraints.append(ppr)
  

# Add all selected pprs to model
for ppr in ppr_restraints:
  ppr.add_to_model()
  restraints.append(ppr)

###########################
# Time everything
###########################
'''
t_cr = timeit.timeit(run_cr,number=1)
t_xl = timeit.timeit(run_xl,number=1)
t_Nic96 = timeit.timeit(run_ppr_Nic96,number=1)
t_Nup53 = timeit.timeit(run_ppr_Nup53,number=1)
t_Nup59 = timeit.timeit(run_ppr_Nup59,number=1)
t_Nup100 = timeit.timeit(run_ppr_Nup100,number=1)
t_Nup116 = timeit.timeit(run_ppr_Nup116,number=1)
t_Nup145N = timeit.timeit(run_ppr_Nup145N,number=1)


print('cr, xl, Nip96, Nup53, Nup59, Nup100, Nup116, Nup145N',t_cr,t_xl,t_Nic96,t_Nup53,t_Nup59,t_Nup100,t_Nup116,t_Nup145N)
print('scores', cr.evaluate(), xl.evaluate(),
      ppr_Nic96.evaluate(), ppr_Nup53.evaluate(), ppr_Nup59.evaluate(),
      ppr_Nup100.evaluate(), ppr_Nup116.evaluate(),ppr_Nup145N.evaluate())
'''
############################
# Scoring function
############################
RS = get_restraint_set(mdl)
sf = IMP.core.RestraintsScoringFunction(RS)
#print('Total score: ', sf.evaluate(False))

############################
# TEST
############################

perms = {('Nup100', 0): [(0,)], ('Nup100', 1): [(1,)]}
prots = list(perms.keys())
iters = list(itertools.product(*perms.values()))[0]
print('prots, iters', prots, list(iters))
start_res = [[560,560],[560,560]]
Tup =  [(p[0], p[1], ii) for p, id in zip(prots, iters) for ii in id]
permutation = list(itertools.chain(*iters))

print('Tup', Tup, permutation)
SEs_ordered = [SEs[i] for i in permutation]

# Finally, evaluate scores
out_perm = open('test.dat','w')
scores_all = eval_threading(Tup, SEs_ordered, start_res, out_perm)
print(scores_all)

print(IMP.atom.show_with_representations(root_hier))
exit()
#zero_coordinates_SEM(SEs_ordered)
############################
# Evaluate all models
############################
for assign, perms in sel_assign.items():
  print('assign', assign, perms)
  permutation = [v[0][0] for k,v in perms.items()]
  perm_label = ''.join([str(p) for p in permutation])
  str_prots = '_'.join([str(a[0])+'.'+str(a[1]) for a in assign])
  out_assign = open(f'enumerate_{str_prots}_{perm_label}.dat','ab')
  # Same prot
  if assign[0][0]==assign[1][0]:
    start_res = [[i,i] for i in range(chains[assign[0][0]][1],chains[assign[0][0]][2],1)]
  else:
    l1 = [i for i in range(chains[assign[0][0]][1],chains[assign[0][0]][2],1)]
    l2 = [i for i in range(chains[assign[1][0]][1],chains[assign[1][0]][2],1)]
    start_res = [[i,j] for i,j in itertools.product(l1,l2)]
  # Get things for scoring
  prots = list(perms.keys())
  iters = list(itertools.product(*perms.values()))[0]
  Tup = [(p[0], p[1], ii) for p,id in zip(prots, iters) for ii in id]
  print('Assign, perms, Tup, start_res length:', assign, perms, Tup, len(start_res))
  SEs_ordered = [SEs[i] for i in permutation]
  # Finally, evaluate scores
  scores_all = eval_threading(Tup, SEs_ordered, start_res, out_assign)
  out_assign.close()
  
exit()






