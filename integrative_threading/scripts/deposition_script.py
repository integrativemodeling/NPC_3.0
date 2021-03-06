# This example demonstrates the use of the Python IHM library to generate
# an mmCIF file for a very simple integrative docking study. Two subunits,
# A and B, each of which is fitted against small angle X-ray (SAXS) data, are
# docked together into a complex, AB, which is fitted against an electron
# microscopy density map.
import sys
sys.path.append("./python-ihm")
import ihm
import ihm.location
import ihm.dataset
import ihm.representation
import ihm.restraint
import ihm.protocol
import ihm.analysis
import ihm.model
import ihm.dumper
import ihm.startmodel
import ihm.cross_linkers
import ihm.reference
import Bio.PDB
import Bio.SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import glob
import os
import itertools
from DisorderedCrossLink import *

sys.path.append('../../utils/')
import make_archive

# The system was represented as a bead model with one residue per bead
class StartingModel(ihm.startmodel.StartingModel):
    def __init__(self, pdb_file, **kwargs):
        super(StartingModel, self).__init__(**kwargs)

        self.pdb_file = pdb_file
        self.asym_unit = kwargs['asym_unit']
        

    def get_atoms(self):
        p = Bio.PDB.PDBParser()
        s = p.get_structure('rep', self.pdb_file)
        for model in s:
            for chain in model.get_chains():
                if chain.get_id()== self.asym_unit.id:
                    print('****', chain.get_id())
                    nasym = asym_mapping[str(chain.get_id())]
                    asym = system.asym_units[nasym]
                    for residue in chain.get_residues():
                        for atom in residue:
                            coord = atom.get_vector()
                            yield ihm.model.Atom(
                                asym_unit=asym, seq_id=residue.get_id()[1],
                                atom_id="CA", type_symbol="C",
                                x=coord[0], y=coord[1], z=coord[2])
                    break

class Model(ihm.model.Model):
    """Pass a BioPython model through to IHM"""
    def __init__(self, file_name, asym_units, **kwargs):
        super(Model, self).__init__(**kwargs)
        self.file_name = file_name
        self.asym_units = asym_units
        self.structured_residues = []

    def get_residues(self):
        # Use BioPython to read the structure from a PDB file, and then yield
        # a set of ihm.model.Residue objects
        p = Bio.PDB.PDBParser()
        s = p.get_structure('rep', self.file_name)

        for model in s:
            for chain in model.get_chains():
                if chain.get_id() in asym_mapping.keys():
                    nasym = asym_mapping[str(chain.get_id())]
                    asym = system.asym_units[nasym]
                    for residue in chain.get_residues():
                        for atom in residue:
                            self.structured_residues.append(residue.get_id()[1])
                            #print(residue, atom)
                            coord = atom.get_vector()
                            self.add_sphere(ihm.model.Sphere(asym_unit=asym, seq_id_range=(residue.get_id()[1], residue.get_id()[1]),
                                                             x=coord[0], y=coord[1],
                                                             z=coord[2], radius=1.0, rmsf=0.0))


    def get_unmodeled_residue_coordinates(self, ps_res):
      # ps_res = Residue to create the pseudo-site for
      #
      #
      # ep1 = Residue of nearest N-terminal endpoint
      # ep2 = Residue of nearest C-terminal endpoint
      # From the threading model, calculate the position of C'v for all unmodeled residues
      # Returns the xyz for each endpoint plus the 
      xyz1 = ep1[0].get_vector()
      xyz2 = ep2[0].get_vector()

      cv_prime = None
      return coords

    def create_pseudo_site(self, residue):
      # Given a residue, find its position according to the model and 
      # return an ihm.restraint.CrossLinkPseudoSite
      coords = [0,0,0]
      ps = ihm.restraint.PseudoSite(coords[0], coords[1], coords[2])
      xlps = ihm.restraint.CrossLinkPseudoSite(ps, self)
      return xlps

    def create_pseudocrosslink(self, expxl, asym1, asym2, distance, resnum1, resnum2):
      # Each crosslink will have two coordinates. These coordinates are either from a modeled residue
      # in which case, the crosslink coordinates are simple the coordiantes of those residues.
      #
      # For those residues that are unmodeled, we must calculate a pseudosite for the residue. This pseudosite
      # will be different for each crosslink and depends on the position of the nearest modeled residues
      xyz1, xyz2 = get_pseudo_site_endpoints(resnum1, resnum2, atoms)
      xlps1 = ihm.restraint.CrossLinkPseudoSite(ihm.restraint.PseudoSite(xyz1[0], xyz1[1], xyz1[2]), self)
      xlps2 = ihm.restraint.CrossLinkPseudoSite(ihm.restraint.PseudoSite(xyz2[0], xyz2[1], xyz2[2]), self)
      rxl = ihm.restraint.ResidueCrossLink(expxl, asym1, asym2, distance, pseudo1=xlps1, pseudo2=xlps2)
      return rxl 

def add_model_cross_links(m, xlrs):
    """Add cross-link information for a given model (the endpoint
       coordinates are model-dependent)."""

    psf = PseudoSiteFinder(m)
    for r in xlrs:
        for xl in r.cross_links:
            ex_xl = xl.experimental_cross_link
            res1 = ex_xl.residue1.seq_id
            res2 = ex_xl.residue2.seq_id
            # Make a function to get this pseudosite
            ps1, ps2, form = psf.get_pseudo_site_endpoints(res1, res2)

            # Add pseudo-site to the existing restraint
            if xl.pseudo1 is None:
                xl.pseudo1 = []
            if xl.pseudo2 is None:
                xl.pseudo2 = []
            xl.pseudo1.append(ps1)
            xl.pseudo2.append(ps2)
            
######################
# MAIN
######################

chain_mapping = {'Nic96':['Q','R','S','T'],
                 'Nup157':['Z','1'],
                 'Nup170':['Y','0'],
                 'Nup53':['U','W'],
                 'Nup59':['V','X'],
                 'Nsp1':['A','D','G','J'],
                 'Nup57':['B','E','H','K'],
                 'Nup49':['C','F','I','L'],
                 'Nup188':['N','P'],
                 'Nup192':['M','O']}
                 #'UNK':['5','6']}


# First, we create a system, which contains everything we know about the
# modeling. A single mmCIF file can contain multiple Systems, but in most
# cases we use just one:

data_dir='../../data/'

title = ("Comprehensive structure and functional adaptations of the yeast nuclear pore complex.")
system = ihm.System(title=title)

system.citations.append(ihm.Citation(
          pmid='34982960', title=title,
          journal="Cell", volume=185, page_range=(361-378),
          year=2022,
          authors=['Akey CW', 'Singh D', 'Ouch C', 'Echeverria I',
                   'Nudelman I', 'Varberg JM', 'Yu Z', 'Fang F',
                   'Shi Y', 'Wang J', 'Salzberg D', 'Song K', 'Xu C',
                   'Gumbart JC', 'Suslov S', 'Unruh J', 'Jaspersen SL',
                   'Chait BT', 'Sali A', 'Fernandez-Martinez J',
                   'Ludtke SJ', 'Villa E', 'Rout MP'],
          doi='10.1016/j.cell.2021.12.015'))


####################
# SOFTWARE
####################

# RaptorX was used to provide a prediction of the secondary structure
system.software.append(ihm.Software(
          name='PSIPRED', classification='secondary structure prediction',
          description='Protein secondary structure prediction based on '
                      'position-specific scoring matrices',
          version='4.0',
          location='http://bioinf.cs.ucl.ac.uk/psipred/'))

# We used various tools from IMP
imp_software = ihm.Software(
          name="Integrative Modeling Platform (IMP)",
          version="2.2",
          classification="integrative model building",
          description="integrative model building",
          location='https://integrativemodeling.org')
system.software.append(imp_software)


####################
# Info about system
####################

# Link to this script
deposition_script = ihm.location.WorkflowFileLocation(
        "deposition_script.py",
        details="Deposition script")
system.locations.append(deposition_script)

# Add uniprot of proteins
lpep = ihm.LPeptideAlphabet()
d = 'Construct used for Cryo-EM'
sd_Nup53 = [ihm.reference.SeqDif(322, lpep['L'], lpep['I'], details=d)]
sd_Nup192 = [ihm.reference.SeqDif(1472, lpep['E'], lpep['A'], details=d)]


Uniprot = {'Nup157':('P40064',[],[]),
           'Nup170':('P38181',[],[]),
           'Nup53':('Q03790',sd_Nup53,[[1,475,1,475]]),
           'Nup59':('Q05166',[],[]),
           'Nsp1':('P14907',[],[]),
           'Nup57':('P48837',[],[]),
           'Nup49':('Q02199',[],[]),
           'Nup188':('P52593',[],[]),
           'Nup192':('P47054',sd_Nup192,[[1,1501,1,1501]]),
           'Nic96':('P34077',[],[[121,838,121,838]])}
           #'Nup100':('Q02629',[],[]),
           #'Nup116':('Q02630',[],[]),
           #'Nup145N':('P49687',[],[[1,605,1,605]])}

# Add all the entities in the model
asym_mapping = {}
k = 0
all_prots = {}
assemblies = []

for prot, (entry, sd, limits) in Uniprot.items():
    for record in Bio.SeqIO.parse(f"{data_dir}/{prot}.fasta", "fasta"):
        if record.name == prot:
            sequence = record.seq
    entity = ihm.Entity(sequence, description=prot)
    system.entities.append(entity)
    chains = chain_mapping[prot]
    for i, c in enumerate(chains):
        asym = ihm.AsymUnit(entity, id=c)
        system.asym_units.append(asym)
        asym_mapping[c] = k
        k += 1
        assembly = ihm.Assembly((asym,), name= f'{prot}.{i}')
        assemblies.append(assembly)
    ref = ihm.reference.UniProtSequence.from_accession(entry)
    for seg in limits:
        ref.alignments.append(ihm.reference.Alignment(
            db_begin=seg[0], db_end=seg[1], entity_begin=seg[2], entity_end=seg[3], seq_dif=sd))
    
    asym.entity.references.append(ref)
        
modeled_assembly = ihm.Assembly(tuple(system.asym_units), name='Modeled assembly')

# Define references sequences
for prot, (entry, sd, limits) in Uniprot.items():
    chains = chain_mapping[prot]
    asym_indexes = [asym_mapping[c] for c in chains]
    ref = ihm.reference.UniProtSequence.from_accession(entry)
    for seg in limits:
        ref.alignments.append(ihm.reference.Alignment(
            db_begin=seg[0], db_end=seg[1], entity_begin=seg[2], entity_end=seg[3], seq_dif=sd))
    for id in asym_indexes[0:1]:
        system.asym_units[id].entity.references.append(ref)
        

# Next, we group asymmetric units (and/or entities) into assemblies.
assembly = ihm.Assembly([asym], name='Threading Model')

####################
# DATA
####################

# XLs data from the 2018 paper
l_XLs = ihm.location.InputFileLocation(
    "../../data/XLs_all_2020_sel_Nic96_Nup53_Nup59_Nup100_Nup116_Nup145N.csv",
    details="DSS Crosslink File")
XLs_dataset = ihm.dataset.CXMSDataset(l_XLs)
XLs_restraint = ihm.restraint.CrossLinkRestraint(XLs_dataset, ihm.cross_linkers.dss)
xlr_dist = 32

all_prots = {e.description:e for e in system.entities}
# Add all experimentally-determined cross-links
with open(XLs_restraint.dataset.location.path) as fh:
    distance = ihm.restraint.UpperBoundDistanceRestraint(xlr_dist)
    for line in fh:
        if line[0]=='#':
            continue
        else:
            vals=line.split(',')
            prot1=vals[0]
            resi1=int(vals[1])
            prot2=vals[2]
            resi2=int(vals[3])
            
            if prot1 in all_prots.keys() and prot2 in all_prots.keys():
                ex_xl = ihm.restraint.ExperimentalCrossLink(
                    all_prots[prot1].residue(resi1),all_prots[prot2].residue(resi2))
                XLs_restraint.experimental_cross_links.append([ex_xl])
                chains1 = chain_mapping[prot1]
                chains2 = chain_mapping[prot2]
                for ch1, ch2 in itertools.product(chains1, chains2):
                    rcl = ihm.restraint.ResidueCrossLink(
                        ex_xl, asym1=system.asym_units[asym_mapping[ch1]],
                        asym2=system.asym_units[asym_mapping[ch2]], distance=distance)
                    XLs_restraint.cross_links.append(rcl)


# EM dataset
l_EM = ihm.location.EMDBLocation('EMDB-24232')
EM_dataset = ihm.dataset.EMDensityDataset(l_EM)
EM_restraint = ihm.restraint.EM3DRestraint(EM_dataset, modeled_assembly)

# PDB data
l_PDB = ihm.location.PDBLocation('7N85')
PDB_dataset = ihm.dataset.PDBDataset(l_PDB)

print(tuple(system.asym_units))

print(type(system.asym_units), system.asym_units)

# Group all together
all_datasets = ihm.dataset.DatasetGroup((XLs_dataset,
                                         EM_dataset,
                                         PDB_dataset))

system.restraints.extend((XLs_restraint,EM_restraint))

####################
# MODELING
####################
protocol = ihm.protocol.Protocol(name='Modeling')
modeling_script = ihm.location.WorkflowFileLocation(
        "mod_SSE_Nic96_copy0_mdff_split.py",
        details="Main modeling script")

protocol.steps.append(ihm.protocol.Step(
                        assembly=modeled_assembly,
                        dataset_group=all_datasets,
                        method='Enumeration',
                        name='Production sampling',
                        num_models_begin=0,
                        num_models_end=1200, multi_scale=False,
                        script_file=modeling_script,
                        software=imp_software))
analysis = ihm.analysis.Analysis()
analysis.steps.append(ihm.analysis.FilterStep(
    feature='energy/score', num_models_begin=1200, num_models_end=5000,
    assembly=modeled_assembly,
    details="Filtering by connectivy satisfaction"))

####################
# Starting model
####################

ST = {}
pdb_ini = '../results/yeast-spoke-thread2-mdff-step13_initial.pdb'
for k,v in asym_mapping.items():
    ST[v]=StartingModel(pdb_ini, asym_unit=system.asym_units[v],dataset=PDB_dataset, asym_id=k)
    ST[v].get_atoms()
    

####################
# Add coordinates
####################

rr = [ihm.representation.ResidueSegment(a, rigid=True, primitive="sphere",starting_model=ST[i]) for i,a in enumerate(system.asym_units)]
rep = ihm.representation.Representation(rr)

# Add ensemble members
models = []
mgs = []
pdbs = glob.glob('../results/yeast-spoke-thread2-mdff-step13_ensemble_*.pdb')
for n, pdb_file in enumerate(pdbs[0:10]):
    print('pdb_file',pdb_file)
    m = Model(assembly = modeled_assembly, protocol=protocol, representation=rep,
              file_name=pdb_file, asym_units=[asym],
              name="Example model %d" % (n))
    m.get_residues()
    add_model_cross_links(m, [XLs_restraint])
    models.append(m)

    r = ihm.location.Repository(
            doi="10.5281/zenodo.5662389",
            url="https://zenodo.org/record/5662389/"
                "files/models_threading.tar.gz")
    
    l_ensemble = ihm.location.OutputFileLocation(path=None, repo=r,
                        details="All models in cluster")

model_group = ihm.model.ModelGroup(models, name="Cluster 1")


me = ihm.model.Ensemble(model_group, len(pdbs))
                        #post_process=protocol.analyses[-1],
                        #file=l_ensemble, name="Cluster %d" % cluster)
system.ensembles.append(me)


# Groups are then placed into states
state = ihm.model.State([model_group])
system.state_groups.append(ihm.model.StateGroup([state]))

# Rewrite local paths to point to Zenodo
repos = []
for subdir, zipname in make_archive.ARCHIVES.items():
    print('subdir', subdir, zipname)
    repos.append(ihm.location.Repository(
        doi="10.5281/zenodo.5662389", root="%s" % subdir,
        url="https://zenodo.org/record/5662389/files/%s.zip" % zipname,
        top_directory=os.path.basename(subdir)))

system.update_locations_in_repositories(repos)

# Once the system is complete, we can write it out to an mmCIF file:
with open('threading_NPC_inner_ring.cif', 'w') as fh:
    ihm.dumper.write(fh, [system])

exit()

