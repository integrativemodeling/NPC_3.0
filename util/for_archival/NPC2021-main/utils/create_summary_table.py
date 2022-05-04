# Load rmf3 file and identify RBs, flexible regions

import IMP
import IMP.rmf
import RMF
import IMP.pmi
import IMP.pmi.tools

import ihm
import ihm.reader

import pandas as pd
import utils
import glob
import os
import numpy as np
from collections import OrderedDict 

class get_input_information(object):
    def __init__(self,
                 mmcif_file):
        
        self.mmcif_file = mmcif_file
        self.datasets = {}
        self.entities = {}

    def get_databases(self):
        model = ihm.model.Model
        system, = ihm.reader.read(open(self.mmcif_file),
                                  model_class=model)

        dbs =  system.orphan_datasets
        for db in system.orphan_datasets:
            print(db)
            self.datasets[db._id] = []
            if db.__class__ == ihm.dataset.ComparativeModelDataset:
                self.datasets[db._id].append(db.__class__)
            elif db.__class__ == ihm.dataset.PDBDataset:
                self.datasets[db._id].append(db.__class__)
                try:
                    self.datasets[db._id].append(db.location.access_code)
                except:
                    print('No access code for:', db._id, db)
            else:
                self.datasets[db._id].append(db.__class__)

        for es in system.entities:
            self.entities[es._id] = [es.description]
                
        mds = system.orphan_starting_models
        for md in system.orphan_starting_models:
            if type(md.dataset) == ihm.dataset.ComparativeModelDataset:
                if len(md.templates)> 0:
                    print('ignacia', md.templates)
                    templates = md.templates
                    for t in templates:
                        db_id = t.dataset._id
                        chain  = t.asym_id
                        seq_range = t.seq_id_range
                        seq_ident = float(t.sequence_identity)
                        # Get pdb access code of template
                        print('aaaaaaa', self.datasets[db_id])
                        if len(self.datasets[db_id]) > 1:
                            pdb_template = self.datasets[db_id][1]
                        else:
                            pdb_template = 'NA'
                        print(md._id, md.dataset.data_type, pdb_template, chain, seq_range, seq_ident)
                        print('entities', self.entities.keys())
                        if md._id in self.entities:
                            self.entities[md._id] += [db_id, md.dataset.data_type, pdb_template, chain, seq_range, seq_ident]
                else:
                    db_id = md.dataset._id
                    if len(self.datasets[db_id]) > 1:
                        pdb_template = self.datasets[db_id][1]
                    else:
                        pdb_template = 'NA'
                    self.entities[md._id] += [db_id, self.datasets[db_id][0].data_type, pdb_template, 'NA', 'NA', 'NA']
            elif type(md.dataset)== ihm.dataset.PDBDataset:
                db_id = md.dataset._id
                self.entities[md._id] += [db_id, md.dataset.data_type,  md.dataset.location.access_code, md.asym_id]
            else:
                db_id = md.dataset._id
                self.entities[md._id] += [db_id, md.dataset, type(md.dataset)]
                
    def get_dictionaries(self):
        # FIX -  read mixed x-ray, comparative
        self.get_databases()
    
        prior_models = []
        for k, p in self.entities.items():
            if len(p)>1:
                if p[2] == 'Comparative model':
                    text = f'{p[0]}: comparative model, template {p[3]}:{p[4]}'

                elif p[2] == 'Experimental model':
                    text = f'{p[0]}: experimental structure, {p[3]}:{p[4]}'
            else:
                text = f'{p[0]}: sequence '

            prior_models.append(text)
            
        temp = OrderedDict()
        temp['Prior models'] = prior_models
        temp['Physical principles and statistical preferences'] = ['Excluded volume', 'Sequence connectivity']

        return temp
    
class get_representation(object):
    def __init__(self,
                 clustering_dir,
                 rmf3=None):
        
        self.clustering_dir = clustering_dir

        print('z')
        if not rmf3:
            rmf3 = os.path.join(self.clustering_dir,'cluster.0/cluster_center_model.rmf')

        self.components = []
        self.struct_components = {}
        self.flex_components = {}
        self.rb_components = {}
        
        self._read_representation(rmf3)
        self._get_rigid_bodies()
        self._unique_components()
        self._get_resolutions(rmf3)
        
    def _read_representation(self,rmf3):
        m = IMP.Model()
        rh_ref = RMF.open_rmf_file_read_only(rmf3)
        h_ref = IMP.rmf.create_hierarchies(rh_ref, m)[0]
        IMP.rmf.load_frame(rh_ref, 0)

        n_struct, n_flex = 0, 0
        for state in h_ref.get_children():
            for component in state.get_children():
                n_residues_prot = 0
                self.components.append(component.get_name())
                for fragment in component.get_children():
                    leaves = IMP.atom.get_leaves(fragment)
                    residues = IMP.atom.Fragment(fragment).get_residue_indexes()
                    residue_range = str(residues[0])+'-'+str(residues[-1])
                    n_res = residues[-1] - residues[0]
                    if IMP.core.RigidMember.get_is_setup(leaves[0]) ==  False:
                        n_flex += n_res
                        if component.get_name() in self.flex_components.keys():
                            self.flex_components[component.get_name()] += ', '+residue_range
                        else:
                            self.flex_components[component.get_name()] = residue_range
                    else:
                        n_struct += n_res
                        for res in fragment.get_children():
                            p = res.get_children()[0]
                            rb = IMP.core.RigidMember(p).get_rigid_body()
                            if rb in self.rb_components.keys():
                                self.rb_components[rb] += [(component.get_name() ,residue_range)]
                            else:
                                self.rb_components[rb] = [(component.get_name() ,residue_range)]

        # Order by name
        self.components = np.sort(self.components)
        
        # Structural coverage
        self.struc_coverage = 100.0* n_struct/(n_struct + n_flex)
        
        import itertools
        for k1, k2 in itertools.combinations(self.rb_components.keys(),2):
            print(k1, k2, k1==k2)

        

    def _get_rigid_bodies(self):
        
        self.RBs = []
        for i, RB in enumerate(self.rb_components.items()):
            all_frags = ''
            for rb_member in RB[1]:
                all_frags += str(rb_member[0]) +'$_{'+rb_member[1]+'}$,'
                if rb_member[0] in self.struct_components.keys():
                    self.struct_components[rb_member[0]] += ', '+rb_member[1]
                else:
                    self.struct_components[rb_member[0]] = rb_member[1]
            self.RBs.append('RB'+str(i+1)+': '+all_frags[0:-1])

    
    def _unique_components(self):

        count_components = {}
        for c in self.components:
            if '.' in c:
                name = c.split('.')[0]
                if name in count_components.keys():
                    count_components[c] += 1
                else:
                    count_components[c] = 1
                
            else:
                count_components[c] = 1

        self.composition = [k +': ' + str(v) for k,v in count_components.items()]
            
    def _get_resolutions(self, rmf3):

        resolution_one = False
        resolutions_flex = []
        resolutions_struct = []
        
        m = IMP.Model()
        rh_ref = RMF.open_rmf_file_read_only(rmf3)
        h_ref = IMP.rmf.create_hierarchies(rh_ref, m)[0]
        IMP.rmf.load_frame(rh_ref, RMF.FrameID(0))

        hiers =  IMP.pmi.tools.input_adaptor(h_ref,
                                             IMP.atom.ALL_RESOLUTIONS,
                                             flatten=True)
        for h in hiers:
            if IMP.core.RigidMember.get_is_setup(h) ==  False:
                r  = len(IMP.atom.Fragment(h).get_residue_indexes())
                if r > 0 and r not in resolutions_flex:
                    resolutions_flex.append(r)
            else:
                r  = len(IMP.atom.Fragment(h).get_residue_indexes())
                if r > 0 and r not in resolutions_struct:
                    resolutions_struct.append(r)
                else:
                    if IMP.atom.Residue(h):
                        resolution_one = True

                        
        if len(resolutions_struct)>0:
            rs = np.max(resolutions_struct)
        else:
            rs = 0
        # FIX    
        if resolution_one and rs > 1:
            self.resolution_structured = '1 [R1], '+ str(rs)+' [R'+str(rs)+'] residues per bead'
        elif resolution_one and rs == 0:
            self.resolution_structured = '1 [R1] residue per bead'
        elif not resolution_one and rs > 1:
            self.resolution_structured = str(rs)+' [R'+str(rs)+'] residues per bead'
        else:
            self.resolution_structured = 'None'


        if len(resolutions_flex):
            rf = np.max(resolutions_flex)
        else:
            rf  = 0
        if rf> 0:
           self.resolution_flexible =  str(rf)+' [R'+str(rf)+'] residues per bead'
        else:
           self.resolution_flexible =  None 
            
    def get_dictionaries(self):
        
        temp = OrderedDict()
    
        temp['Composition (number of copies)'] = self.composition
        temp['Atomic (structured) components'] = [k+': '+self.struct_components[k] for k in sorted(self.struct_components)]
        if len(self.flex_components.items())> 0:
            temp['Unstructured components'] = [k+': '+self.flex_components[k] for k in sorted(self.flex_components)]
        else:
            temp['Unstructured components'] = ['None']
        temp['Resolution of structured components'] = [self.resolution_structured]
        temp['Resolution of unstructured components'] = [self.resolution_flexible] 
        temp['Structural coverage'] = [str(round(self.struc_coverage,2))+' \%']
        temp['Rigid body (RB) definitions'] = [rb for rb in self.RBs]

        #FIX
        temp['Spatial restraints encoded into scoring function'] = ['Excluded volume; applied to the R1 representation',
                                                                    'Sequence connectivity; applied to the R1 representation']
    
       
        return temp
        
class read_modeling_information(object):
    def __init__(self,
                 modeling_script,
                 analysis_directory,
                 clustering_directory):

        self.modeling_script = modeling_script
        self.clustering_dir = clustering_directory
        self.analysis_directory = analysis_directory

        self.clusters = []
        
        S = pd.read_csv(os.path.join(self.analysis_directory,'summary_sampling_information.csv'),
                        sep=',', index_col=0, header=None)
        self.S = S.T
        self.parse_modeling_script()

        self.get_KS_stats()

        self.get_sampling_precision()
        self.get_clusters_precision()
        self.get_cluster_population()

        self._compute_cross_correlation()
        
    def parse_modeling_script(self):
        self.restraints = []
        self.max_rb_trans = 4.0
        self.max_rb_rot = 0.4
        self.max_bead_trans = 4.0
        for line in open(self.modeling_script):
            if 'IMP.pmi.restraints.' in line and 'import' not in line:
                r = line.split('IMP.pmi.restraints.')[1].split('(')[0]
                if '.' in r:
                    r = r.split('.')[-1]
                self.restraints.append(r)
            if 'replica_exchange_minimum_temperature' in line:
                self.min_temp = line.split('replica_exchange_minimum_temperature')[1].split('=')[1].split(',')[0]
            else:
                self.min_temp = 1.0

            if 'replica_exchange_maximum_temperature' in line:
                self.max_temp = line.split('replica_exchange_maximum_temperature')[1].split('=')[1].split(',')[0]
            else:
                self.max_temp = 2.5
                
            if 'max_rb_trans' in line:
                self.max_rb_trans = line.split('max_rb_trans')[1].split('=')[1].split(',')[0].replace(')','')
                if type(self.max_rb_trans) is str:
                    self.max_rb_trans = 4.0
                
            if 'max_rb_rot' in line:
                self.max_rb_rot = line.split('max_rb_rot')[1].split('=')[1].split(',')[0].replace(')','')
                if type(self.max_rb_rot) is str:
                    self.max_rb_rot = 0.4
            if 'max_bead_trans' in line:
                self.max_bead_trans = line.split('max_bead_trans')[1].split('=')[1].split(',')[0].replace(')','')
                if type(self.max_bead_trans) is str:
                    self.max_bead_trans = 4.0
                    
                
        print(self.restraints)
        
        
    def get_dictionaries_sampling(self):

        temp = OrderedDict()
        temp['Sampling method'] = ['Replica Exchange Gibbs sampling, based on Metropolis Monte Carlo']
        temp['Replica exchange temperature range'] = [f'{self.min_temp} - {self.max_temp}']
        temp['Number of replicas'] = self.S['Number_of_replicas'].values[0]
        temp['Number of runs'] = self.S['Number_of_runs'].values[0]
        temp['Number of structures generated'] = self.S['N_total'].values[0]
        #temp['Movers for rigid bodies'] = [f'Random translation up to {self.max_rb_trans} \AA', f'Random rotation up to {self.max_rb_rot} radians']
        temp['Movers for flexible string of bead'] = [f'Random translation up to {self.max_bead_trans} \AA']

        return temp

    def get_dictionaries_models(self):

        temp = OrderedDict()
        temp['Number of models after equilibration'] = self.S['N_total'].values[0]
        temp['Number of models that satisfy the input information'] = int(float(self.S['N_selected'].values[0]))
        temp['Number of structures in samples A/B'] = [str(self.S['N_sample_A'].values[0])+'/'+ str(self.S['N_sample_B'].values[0])]
        temp['p-value of non-parametric Kolmogorov-Smirnov two-sample test'] = [str(round(float(self.KS_0),3))+' (threshold p-value $>$ 0.05)']
        temp['Kolmogorov-Smirnov two-sample test statistic, D'] = [str(round(float(self.KS_1),2))]
        
 
        return temp

    def get_dictionaries_clustering(self):

        # Report only clusters with 2% or more models
        self.cluster_PRE = self.cluster_PRE[0:len(self.perc_POP)]
        
        temp = {}
        
        temp['Sampling precision'] = [str(round(float(self.sampling_precision),2))+' \AA ']
        temp["Homogeneity of proportions $\chi^2$ test (p-value)/Cramer’s V value"] = [f"{self.p_value}/{self.cramers_V} (thresholds: p-value$>$0.05 OR Cramer's V$<$0.1)"]
        temp['Number of clusters'] = [str(len(self.perc_POP))]
        temp['Cluster populations'] =  self.perc_POP
        temp['Cluster precisions'] =  self.cluster_PRE
        temp['Average cross-correlation between localization probability densities of samples A and B'] = ['cluster '+ str(n+1)+': '+ str(cluster['average_cc'])  for n, cluster in enumerate(self.clusters[0:len(self.perc_POP)])]
        
        return temp
    
    def get_KS_stats(self):
        file_in = glob.glob(self.clustering_dir+'/*'+'KS_Test.txt')
        for line in open(file_in[0]):
            vals = line.split()
            self.KS_0 = vals[0]
            self.KS_1 = vals[1]    
        
    def get_sampling_precision(self):
        file_in = glob.glob(self.clustering_dir+'/*'+'Sampling_Precision_Stats.txt')
        print(file_in)
        with open(file_in[0], 'r') as f:
            for line in f:
                if line.startswith("Sampling "):
                    info = next(f)
                    break
        info = info.split()
        self.sampling_precision = info[0]
        self.p_value = info[1]
        self.cramers_V = info[2]
        self.percent = info[3]

        print('Sampling precision', self.sampling_precision, self.p_value, self.cramers_V, self.percent)

    def get_clusters_precision(self):
        file_in = glob.glob(self.clustering_dir+'/*'+'Cluster_Precision.txt')
        c = 0
        for line in open(file_in[0],'r'):
            vals = line.split()
            self.clusters.append({})
            self.clusters[c]['precision'] = vals[11]
            c += 1

        self.cluster_PRE = ['cluster '+ str(i+1)+' : '+ str(round(float(p['precision']),2)) + ' \AA'
                            for i, p in enumerate(self.clusters)]
        
    def get_cluster_population(self):
        self.total_clustered, c = 0, 0
        file_in = glob.glob(self.clustering_dir+'/*'+'Cluster_Population.txt')
        for line in open(file_in[0],'r'):
            vals = line.split()
            self.clusters[c]['population'] = float(vals[2])
            self.total_clustered +=  float(vals[2])
            c += 1 
            
        # Compute percent
        self.perc_POP = []
        for i, p in enumerate(self.clusters):
            pop = 100*float(p['population'])/float(self.total_clustered)
            if pop > 5.0:
                self.perc_POP.append('cluster '+ str(int(i+1))+' : '+str(round(pop,1))+' $\%$')
        
    def _get_cluster_densities(self):

        for ncluster, cluster in enumerate(self.clusters):
            dens_A = glob.glob(os.path.join(
                self.clustering_dir,'cluster.'+str(ncluster),'Sample_A/*.mrc'))
            dens_B = glob.glob(os.path.join(
                self.clustering_dir,'cluster.'+str(ncluster),'Sample_B/*.mrc'))

            
            
            dens_A = [d for d in dens_A if 'all' not in d.lower()]
            dens_B = [d for d in dens_B if 'all' not in d.lower()]

            print('A, B', dens_A, dens_B)
            
            cluster['densities']  = zip(dens_A, dens_B)
            
        
    def _compute_cross_correlation(self):

        self._get_cluster_densities()

        mrc_rw = IMP.em.MRCReaderWriter()
        ccc = IMP.em.CoarseCC()

        for ncluster, cluster in enumerate(self.clusters):
            cc_vals = []
            for n, (dens_A, dens_B) in enumerate(cluster['densities']):
                dmapA = IMP.em.read_map(dens_A, mrc_rw)
                dmapB = IMP.em.read_map(dens_B, mrc_rw)    
                mass = IMP.atom.get_mass_from_number_of_residues(1000.)
                threshold = IMP.em.get_threshold_for_approximate_mass(dmapA, mass)
                score = ccc.cross_correlation_coefficient(dmapA, dmapB, threshold, True)
                cc_vals.append(score)
            cluster['average_cc'] = round(np.mean(cc_vals), 2)

        print(self.clusters)
            
    def update_mmcif_file(self, mmcif_file):

        model = ihm.model.Model
        system, = ihm.reader.read(open(mmcif_file,'r'), model_class=model)

        protocol = system.orphan_protocols[-1]
        protocol.steps[-1].num_models_end = self.S['N_total']
        
        # Next, filtered  models
        analysis = ihm.analysis.Analysis()
        protocol.analyses.append(analysis)
        analysis.steps.append(ihm.analysis.FilterStep(
            feature='energy/score',
            num_models_begin=self.S['N_total'], num_models_end=self.S['N_selected']))
        

        #tmpd = tempfile.mkdtemp()
        for ncluster, cluster in enumerate(self.clusters):
            e = ihm.model.Ensemble(model_group='m1',
                           #file ,
                           num_models=cluster['population'],
                           post_process=None,
                           name='cluster%s'%(ncluster+1),
                           clustering_method='Hierarchical',
                           clustering_feature='RMSD',
                           precision=cluster['precision'])

            d = ihm.model.LocalizationDensity(
                os.path.join(self.clustering_dir,'cluster.%s/LPD_All.mrc' %ncluster),
                                              asym_unit='cluster%s'%(ncluster + 1))

            e.densities.append(dens)
            system.ensembles.append(e)

            
            #r = ihm.location.Repository(doi="10.5281/zenodo.1445841",
            #                           url="https://zenodo.org/record/1445841/files/cluster%d.dcd"% ncluster)
            #f = ihm.location.OutputFileLocation(path='.',
            #                                    #repo=r,
            #                                    details="All ensemble structures for cluster %d" % ncluster)
            #e = po._add_simple_ensemble(analysis.steps[-1],
            #                            name="Cluster %d" % ncluster,
            #                            num_models=cluster['size'],
            #                            rmsd=cluster['precision'],
            #                            num_models_deposited=1,
            #                            localization_densities={},
            #                            ensemble_file=f)
            # Add localization density 
            #loc = ihm.location.OutputFileLocation(
            #    self.clustering_dir+'/cluster.%s/LPD_All.mrc' % (ncluster))
        #    den = ihm.model.LocalizationDensity(file=loc,
        #                                        asym_unit=po.asym_units['ecm29.0'])
        #    e.densities.append(den)

            # Add one output model
            #rmf_file = fix_rmf_file(self.clustering_dir+'/cluster.%s/cluster_center_model.rmf3' % (ncluster)),
         #                           moldict, tmpd)
            #rh = RMF.open_rmf_file_read_only(rmf_file)
            #IMP.rmf.link_hierarchies(rh, [representation])
            #IMP.rmf.load_frame(rh, RMF.FrameID(0))
            #del rh

        #model = po.add_model(e.model_group)
        #shutil.rmtree(tmpd)
        

class read_validation_information(object):
    def __init__(self,
                 clustering_dir):
        
        self.clustering_dir = clustering_dir

    def get_dictionaries(self):
        temp = OrderedDict()
        temp['Percent of sequence connectivity restraints satisfied per structure'] = ['99 \%']
        return temp

class read_benchmark_information(object):
    def __init__(self,
                 clustering_dir):
        self.clustering_dir = clustering_dir

    def get_dictionaries(self):
        temp = OrderedDict()
        temp['Structural accuracy (95 \% CI)'] = []
        
        return temp
    

class get_software_information(object):
    def __init__(self,
                 mmcif_file):
        self.mmcif_file = mmcif_file

        model = ihm.model.Model
        sys, = ihm.reader.read(open(mmcif_file,'r'), model_class=model)
        all_s = sys.software
        self.software = []
        for s in all_s:
            self.software.append(s.name+', version '+s.version)
    
    def get_dictionaries(self):
        temp = OrderedDict()
        temp['Modeling programs'] = self.software
        
        return temp
    
class get_data_availability(object):
    def __init__(self, mmcif_file):
        self.mmcif_file = mmcif_file
        
    def get_dictionaries(self):
        temp = OrderedDict()
        temp['PDB-dev accesion code'] = ['TBD']
    
        return temp

def dict_to_list(d):
    L = []
    for k,v in d.items():
        if isinstance(v, list):
            L.append([k,v])
        else:
            L.append([k,[v]])
    return L
'''
###########################
# Read rmf3
###########################
modeling_script = '/Users/iecheverria/Dropbox/UCSF/pE-MAP/modeling_pmi2/repo/modeling_rnap/modeling/mod_pemap_rnap.py'
mmcif_file = '/Users/iecheverria/Dropbox/UCSF/pE-MAP/modeling_pmi2/repo/modeling_rnap/RNAP_pEMAP_initial.cif'

analysis_dir = '/Users/iecheverria/Dropbox/UCSF/pE-MAP/modeling_pmi2/results/results_RNAP_XLs/analys/'
clustering_dir = os.path.join(analysis_dir,'clustering_cl5')
rmf3 = os.path.join(clustering_dir,'cluster.0','cluster_center_model.rmf3')

I = get_input_information(mmcif_file)
input_information = I.get_dictionaries()

R = get_representation(clustering_dir)
representation = R.get_dictionaries()

S = read_modeling_information(modeling_script,
                              analysis_dir,
                              clustering_dir)

sampling = S.get_dictionaries_sampling()

samples = S.get_dictionaries_models()

clustering = S.get_dictionaries_clustering()

#S.update_mmcif_file(mmcif_file)

V = read_validation_information(clustering_dir)
validation = V.get_dictionaries()

V = read_benchmark_information(clustering_dir)
benchmark = V.get_dictionaries()

SS = get_software_information(mmcif_file)
software = SS.get_dictionaries()

print(software)
        
D = get_data_availability(clustering_dir)
data_availability = D.get_dictionaries()

################################################
# Edit dictionaries
# Entries is dictionaries can be edited to add
# other custom information
################################################
input_information['Experimental data'] = ['123 pE-MAP derived distance restraints']

representation['Spatial restraints encoded into scoring function'].append('pE-MAP MIC pair-restraints; applied to the R1 representation')

representation['Rigid body (RB) definitions'] = ['RB1: rpb1$_{1-1105}$', 'RB2: rpb1$_{1113-1404}$', 'RB3: rpb2$_{1-1099}$']

sampling['CPU time'] = ['12 hours on 20 processors']

software['Modeling scripts'] = ['https://github.com/salilab/pemap']
software['Homology detection and structure prediction'] = ['HHPred, version 2.0.16']
software['Visualization and plotting'] = ['UCSF Chimera, version 1.10', 'Matplotlib, version 3.0.3 ']

################################################
# Convert ordered dictionaries 
# into lists
################################################
input_information_list = dict_to_list(input_information)
representation_list = dict_to_list(representation)
sampling_list = dict_to_list(sampling)
samples_list = dict_to_list(samples)
clustering_list = dict_to_list(clustering)
validation_list = dict_to_list(validation)
benchmark_list = dict_to_list(benchmark)
software_list = dict_to_list(software)
data_availability_list = dict_to_list(data_availability)


print(sampling_list)


################################################
# Compile all information
# 
################################################
variable_dict = {'complex': 'RNAPII rpb1-rpb2 dimer',
                 'input_information': input_information_list, 
                 'representation': representation_list,
                 'sampling': sampling_list,
                 'samples': samples_list,
                 'clustering':clustering_list,
                 'validation':validation_list,
                 'benchmark':benchmark_list,
                 'software':software_list,
                 'data':data_availability_list}


################################################
# Generate tex, pdf file
# 
################################################
template = utils.get_template('./SI_template.tex')
utils.compile_pdf_from_template(template, variable_dict, './table_SI_RNAP.pdf')

exit()
'''
