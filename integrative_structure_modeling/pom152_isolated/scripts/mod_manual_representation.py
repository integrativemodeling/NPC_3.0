import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.mmcif
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em 
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.npc_restraints
import IMP.bayesianem
import IMP.bayesianem.restraint
import math
import time
import numpy as np
import sys


from sys import exit

##############################
# Restraint to include
##############################

top_dir = '../'

include_XLs = True
include_EM = True
include_intra_XLs = True
include_IG_distance = True
include_Memb_binding = True

w_xls = 2.0
w_em = 50.0
w_msl = 10.0
w_intra_xls = 10.0
w_ig_distance  = 10.0 

# Create System and State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

# Read sequences and create Molecules
seqs = IMP.pmi.topology.Sequences(top_dir+'data/yPom152.fasta')
mol = st.create_molecule("Pom152",sequence=seqs["yPom152"],chain_id='A')

# Add structure. This function returns a list of the residues that now have structure
a1 = mol.add_structure(top_dir+'data/pom152_copies_Ig0_0.pdb',
                          chain_id='A')
a2 = mol.add_structure(top_dir+'data/pom152_E.pdb',
                          chain_id='E')

# Add structured part representation and then build
mol.add_representation(mol.get_atomic_residues(),
                       density_prefix='Pom152',
                       density_residues_per_component=60,
                       density_voxel_size=3.0,
                       resolutions=[1])
mol.add_representation(mol.get_non_atomic_residues(),
                       #density_prefix='Pom152',
                       #density_residues_per_component=20,
                       #density_voxel_size=3.0,
                       resolutions=[10])

print('building clones')

mols = [mol]
chains='BCDEFGHI'
for nc in range(7):
    clone = mol.create_clone(chains[nc])
    mols.append(clone)

hier = s.build()

print(mols)
for m in mols:
    print(m.get_name())

print([ss.get_molecules() for ss in s.get_states()])
    


##############################
# Generate mmcif file
##############################
    
if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(open('Pom152_ring_isolated.cif', 'w'))
    po.system.title = ('Integrative structure determination of Pom152 rings of isolated NPCs')
    s.add_protocol_output(po)

# Create a symmetry constraint
#  A constrant is invariant: IMP will automatically move all clones to match the reference
#  If instead you want some more flexiblity, consider IMP.pmi.restraints.stereochemistry.SymmetryRestraint

##############################
# DOF and symmetry
##############################

dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

domains = [ [260,362],
            [379,472],
            [520,611],
            [616,714],
            [722,818],
            [824,918],
            [931,1026],
            [1036,1141],
            [1150,1229],
            [1244,1337]]

for i, mol in enumerate(mols):
    dof.create_flexible_beads(mol)
    for j, d in enumerate(domains):
        sel = IMP.atom.Selection(hier,
                                 molecule='Pom152',
                                 residue_indexes=range(d[0],d[1]+1),
                                 copy_index=i).get_selected_particles()
    
        rb_movers,rb = dof.create_rigid_body(sel,
                                             name = f'test_{i}_flex')
       
    
# Symmetry

center = IMP.algebra.Vector3D([0,0,0])

# Axial symmetry
rot = IMP.algebra.get_rotation_about_axis([1,0,0],math.pi)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[0],mols[1],transform)

# Radial symmetry
rot = IMP.algebra.get_rotation_about_axis([0,0,1],2*math.pi*(1)/8)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[0],mols[2],transform)
dof.constrain_symmetry(mols[1],mols[3],transform)

rot = IMP.algebra.get_rotation_about_axis([0,0,1],2*math.pi*(2)/8)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[0],mols[4],transform)


rot = IMP.algebra.get_rotation_about_axis([0,0,1],2*math.pi*(6)/8)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[1],mols[5],transform)

rot = IMP.algebra.get_rotation_about_axis([0,0,1],2*math.pi*(7)/8)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[0],mols[6],transform)
dof.constrain_symmetry(mols[1],mols[7],transform)

mdl.update() # propagates coordinates

out = IMP.pmi.output.Output()
out.init_rmf("sym_ini.rmf3", [hier])
out.write_rmf("sym_ini.rmf3")
out.close_rmf("sym_ini.rmf3")

# Some restraints

#############################
# Connectivity restraint
#############################

crs = []
output_objects = []
mols_sim = []
mols_sym = []
for mol in mols:
    copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
    print('copy', copy_n)
    if copy_n<1:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.set_label(mol.get_name()+'.'+str(copy_n))
        cr.add_to_model()
        output_objects.append(cr)
        mols_sim.append(mol)
    else:
        mols_sym.append(mol)

##########################
# EM restraint
##########################
if include_EM:
    densities = IMP.atom.Selection(hier,
                                   representation_type=IMP.atom.DENSITIES).get_selected_particles()

    print(densities, len(densities))
    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                              top_dir+'data/em_data/run8_locres-box-480-pom152-cutout-lp25A-centered_3spokes_dsfact2_ng80.txt',
                                                              scale_target_to_mass=True)
    gem.set_label("EM")
    gem.add_to_model()
    gem.set_weight(w_em)
    output_objects.append(gem)

    #gem.center_model_on_target_density(st)
    t0 = gem.evaluate()
    
    print('Eval. EM at t0: ', t0)


##########################
# Inter-molecular XLs
##########################
resi_XLs = [62, 301, 351]

if include_intra_XLs:
    for r in resi_XLs:
        dist_min = 3.0
        dist_max = 20.0
        ixl = IMP.pmi.restraints.basic.DistanceRestraint(root_hier = hier,
                                                         tuple_selection1=(r,r,'Pom152',0),
                                                         tuple_selection2=(r,r,'Pom152',1),
                                                         distancemin=dist_min,
                                                         distancemax=dist_max,
                                                         label=f'XLs_inter_{r}')
        ixl.set_weight(w_intra_xls)
        ixl.add_to_model()
        output_objects.append(ixl)
        print('Intra molecular XLs:', ixl.get_output())
    
                
##########################
# Distance between IG
# Restraint
##########################

if include_IG_distance:
    IG_list = [[362, 379], [471,520], [611,616], [713,722], [816,824], [916,933], [1025,1038], [1141,1150], [1229,1244] ]
    for rs in IG_list:
        dist_min = 3.0
        dist_max = 20.0
        dr = IMP.pmi.restraints.basic.DistanceRestraint(root_hier = hier,
                                                        tuple_selection1=(rs[0],rs[0],'Pom152',0),
                                                        tuple_selection2=(rs[1],rs[1],'Pom152',0),
                                                        distancemin=dist_min,
                                                        distancemax=dist_max)
        dr.set_weight(w_ig_distance)
        dr.set_label(f'IG_distance_{rs[0]}')
        dr.add_to_model()
        output_objects.append(dr)
        print('IG distance restraint', dr.get_output())

##########################
# Membrane binding
##########################
tor_th      = 45.0
tor_th_ALPS = 12.0
tor_R       = 390.0 + 150.0
tor_r       = 150.0 - tor_th/2.0
tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
msl_sigma   = 1.0
        
###########################
# Membrane binding
###########################

if include_Memb_binding:

    membrane_surface_sel = [(111,194,'Pom152',0)]

    for sel in membrane_surface_sel:
        msl = IMP.pmi.restraints.npc_restraints.MembraneSurfaceLocationRestraint(hier=hier,
                                                                                   protein=sel,
                                                                                   tor_R=tor_R,
                                                                                   tor_r=tor_r_ALPS,
                                                                                   tor_th=tor_th_ALPS,
                                                                                   sigma=msl_sigma,
                                                                                   resolution = 1)
        msl.set_weight(w_msl)
        msl.add_to_model()
        msl.set_label('%s.%s'%(sel[2],sel[0]))
        output_objects.append(msl)
        print('Membrane binding restraint:', msl.evaluate())

###########################
# Randomize configurations
###########################

sel = IMP.atom.Selection(hier,
                         molecule='Pom152').get_selected_particles()

IMP.pmi.tools.shuffle_configuration(sel,
                                    bounding_box=((400, -300, 0), (550, 0, 80)),
                                    avoidcollision_rb=False)
        
############################
# Sampling
############################

num_frames = 100
if '--mmcif' in sys.argv or '--test' in sys.argv:
    num_frames=5

mc1 = IMP.pmi.macros.ReplicaExchange0(mdl,
                                      root_hier=hier,                       
                                      #crosslink_restraints=rmf_restraints,       
                                      monte_carlo_sample_objects=dof.get_movers(),  
                                      global_output_directory="pre_output/",
                                      output_objects=output_objects,
                                      replica_exchange_maximum_temperature=4.0,
                                      monte_carlo_steps=10,
                                      number_of_frames=num_frames,
                                      number_of_best_scoring_models=0)

mc1.execute_macro()
rex1 = mc1.get_replica_exchange_object()

############################
# Excluded Volume
############################

print(mols_sim, mols_sym)

# Create excluded volume for all particles
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols_sim)
evr1.add_to_model()
evr1.set_weight(1.0)
evr1.set_label('intra')
output_objects.append(evr1)

evr2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols_sim,
                                                               other_objects=mols_sym)
evr2.add_to_model()
evr2.set_weight(1.0)
evr2.set_label('inter')
output_objects.append(evr2)

print('Excluded volume:', evr1.evaluate(), evr2.evaluate())

###########################
# Chemical crosslinks
###########################
# INITIALIZE DB
rmf_restraints = []
if include_XLs:
    cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    cldbkc.set_protein1_key("Protein1")
    cldbkc.set_protein2_key("Protein2")
    cldbkc.set_residue1_key("Residue1")
    cldbkc.set_residue2_key("Residue2")
    #cldbkc.set_unique_id_key("UniqueID")
    #cldbkc.set_psi_key("Psi")

    # XLs RESTRAINT
    cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb.create_set_from_file(top_dir+"data/XLs_all_2020_Pom152.csv")

    xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                                database=cldb,
                                                                                resolution=1,
                                                                                length=21.0,
                                                                                slope=0.01)
    xl1.add_to_model()
    xl1.set_weight(w_xls)
    rmf_restraints.append(xl1)
    output_objects.append(xl1)

############################
# Sampling
############################

num_frames = 50000
if '--mmcif' in sys.argv or '--test' in sys.argv:
    num_frames=5


mc2 = IMP.pmi.macros.ReplicaExchange0(mdl,
                                      root_hier=hier,                       
                                      crosslink_restraints=rmf_restraints,       
                                      monte_carlo_sample_objects=dof.get_movers(),  
                                      global_output_directory="output/",
                                      output_objects=output_objects,
                                      replica_exchange_maximum_temperature=3.0,
                                      monte_carlo_steps=10,
                                      number_of_frames = num_frames,
                                      number_of_best_scoring_models=0,
                                      replica_exchange_object = rex1)

mc2.execute_macro()


##############################
# Generate mmcif
##############################  
if '--mmcif' in sys.argv:
    import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    
    fname = top_dir+"data/XLs_all_2020_Pom152.csv"

    # Add publication
    po.system.title = "Integrative structure determination of the NPC"
    #po.system.citations.append(ihm.Citation.from_pubmed_id(32719457))
    
    s = po.system
    print("restraint datasets:", [r.dataset for r in s.restraints])
    # Datasets for XL-MS restraint
    for r in s.restraints:
        if isinstance(r, ihm.restraint.CrossLinkRestraint):
            r.linker = ihm.cross_linkers.dsso
            print("XL-MS dataset at:", r.dataset.location.path)
            print("Details:", r.dataset.location.details)
    
    # Correct number of output models to account for multiple runs
    protocol = s.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 2007800

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 9999 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=200000,
                            num_models_end=9999))
    

    # Create an ensemble for the cluster
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0", num_models=9999,
                                drmsd=8.3, num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    #rh = RMF.open_rmf_file_read_only('../results/clustering_rigid/cluster.0/cluster_center_model.rmf3')
    #IMP.rmf.link_hierarchies(rh, [hier])
    #IMP.rmf.load_frame(rh, RMF.FrameID(0))
    #del rh
    #model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'../results/clustering/cluster.0/LPD_{name}.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
    
    # Add uniprot of proteins
    # name : (uniprot id, mutations, [[db_begin, db_end, entity_begin, entity_end]]
    Uniprot={'Pom152': ('P39685',[],[])}

    #for prot, (entry, sd, limits) in Uniprot.items():
    #    print(prot, entry, sd, limits)
    #    ref = ihm.reference.UniProtSequence.from_accession(entry)
    #    for seg in limits:
    #        ref.alignments.append(ihm.reference.Alignment(
    #            db_begin=seg[0], db_end=seg[1], entity_begin=seg[2], entity_end=seg[3], seq_dif=sd))
            
    #    po.asym_units[prot].entity.references.append(ref)

    # Point to the raw mass spec data and peaklists used to derive the crosslinks.
    #l = ihm.location.PRIDELocation('PXD019338',
    #                               details='All raw mass spectrometry files and '
    #                               'peaklists used in the study')
    #xl1.dataset.add_primary(ihm.dataset.MassSpecDataset(location=l))
        

    # Replace local links with DOIs
    repos = []
    #for subdir, zipname in make_archive.ARCHIVES.items():
    #    print('subdir', subdir)
    #    repos.append(ihm.location.Repository(
    #        doi="10.5281/zenodo.3836213", root="../%s" % subdir,
    #        url="https://zenodo.org/record/3836213/files/%s.zip" % zipname,
    #        top_directory=os.path.basename(subdir)))
    
    po.system.update_locations_in_repositories(repos)
    

    po.flush()


