###################################
# Script to summarize all modeling
# information into one table
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
###################################

import pandas as pd
import glob
import os
import sys
import numpy as np

sys.path.append('../utils')
from create_summary_table import *
import utils

###########################
# Read files
###########################
modeling_script = '../integrative_structure_modeling/pom152_isolated/scripts/mod_manual_representation.py'
mmcif_file = '../integrative_structure_modeling/pom152_isolated/scripts/Pom152_ring_isolated.cif'

analysis_dir = '../integrative_structure_modeling/pom152_isolated/results/'
clustering_dir = os.path.join(analysis_dir,'clustering')
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

#V = read_benchmark_information(clustering_dir)
#benchmark = V.get_dictionaries()

SS = get_software_information(mmcif_file)
software = SS.get_dictionaries()

D = get_data_availability(clustering_dir)
data_availability = D.get_dictionaries()

################################################
# Edit dictionaries
# Entries is dictionaries can be edited to add
# other custom information
################################################
input_information['Experimental data'] = ['1425 DSS and EDC chemical crosslinks']
input_information['Experimental data'].append('Cryo-EM map; EMDB TBD')
input_information['Experimental data'].append('FIB-milled CET map; EMDB TBD')
input_information['Experimental data'].append('$in\,  vitro$ protein-protein interactions')
input_information['Prior models'].append('Structure of the isolated and $in\, situ$ NPCs')
input_information['Physical principles and statistical preferences'].append('Predicted secondary structure (RaptorXProperties)')
input_information['Physical principles and statistical preferences'].append('Predicted transmembrane domains (http://yeastgenome.org and HeliQuest)')

print(representation.keys())


representation['Number of structural elements (SEs)'] = 30
representation['Structural coverage'] = '71 \%'

representation['Resolution of structured components'] = '1 [R1] residue per bead'

representation['Composition (number of copies)'] = 8
representation['Composition (number of copies of Pom152)'] = representation.pop('Composition (number of copies)')

representation['Atomic (structured) components'] = '260-362, 379-472, 520-611, 616-714, 722-818, 824-918, 931-1026, 1036-1141, 1150-1229, 1244-1337'
representation['Atomic (structured) components (Pom152)'] = representation.pop('Atomic (structured) components')

representation['Unstructured components'] = '1-259, 361-378, 471-519, 610-615, 713-721, 819-823, 919-930, 1027-1035, 1142-1149, 1230-1243'
representation['Unstructured components (Pom152)'] = representation.pop('Unstructured components')


representation['Spatial restraints encoded into scoring function (Pom152)'] = representation.pop('Spatial restraints encoded into scoring function')
representation['Spatial restraints encoded into threading scoring function'] = []

representation['Spatial restraints encoded into threading scoring function'].append('Secondary structure restraint; applied to the SEs')
representation['Spatial restraints encoded into threading scoring function'].append('Loop end-to-end disntance restraint; applied to the SEs')
representation['Spatial restraints encoded into threading scoring function'].append('Cross-link restraints; applied to the SEs and R1 representation')

representation['Spatial restraints encoded into scoring function (Pom152)'].append('Cross-link restraints; applied to the R1 representation')
representation['Spatial restraints encoded into scoring function (Pom152)'].append('EM density restraint using Gaussian Mixture Model (GMMs) representations')
representation['Spatial restraints encoded into scoring function (Pom152)'].append('Transmembrane domain restraint (Pom152 (111-200))')
representation['Spatial restraints encoded into scoring function (Pom152)'].append('Peri-nuclear restraint (Pom152 (201-1337))')
representation['Spatial restraints encoded into scoring function (Pom152)'].append('Sequence connectivity and excluded volume')

enumeration = {}
enumeration['Sequences used for threading'] = 'Nic96 (residues 1-205), Nup53 (residues 1-247), Nup59 (residues 1-265), Nup100 (residues 551-815), Nup116 (residues 751-965), and Nup145N (residues 201-458)'

sampling['CPU time'] = ['14 hours on 42 processors']

validation_threading = {'Protein, sequence range, and start residue standard deviation' : ['SE1: Nic96, 15-30 (1)',
                                                                                           'SE2: Nic96, 38-43 (2)',
                                                                                           'SE3: Nic96, 46-61 (0.3)',
                                                                                           'SE4: Nic96, 66-78 (2)',
                                                                                           'SE5: Nic96, 84-99 (1)',
                                                                                           'SE6: Nic96, 118-141 and 147-166 (1.4)',
                                                                                           'SE7: Nic96, 118-153 and 155-165 (0)',
                                                                                           'SE8: undefined']}

validation['Percent of sequence connectivity restraints satisfied per structure'] = ['99/99 \%']
validation['Percent cross-link restraints satisfied by ensemble'] = ['89/89 \%']
validation['Percent of transmembrane domain restraints satisfied by ensemble'] = ['98/99 \%']
validation['Percent of excluded volume restraints satisfied per structure'] = ['99/99 \%']


validation['Percent of sequence connectivity restraints satisfied per structure (isolated/in situ)'] = validation.pop('Percent of sequence connectivity restraints satisfied per structure')
validation['Percent cross-link restraints satisfied by ensemble (isolated/in situ)'] = validation.pop('Percent cross-link restraints satisfied by ensemble')
validation['Percent of transmembrane domain restraints satisfied by ensemble (isolated/in situ})'] =  validation.pop('Percent of transmembrane domain restraints satisfied by ensemble')
validation['Percent of excluded volume restraints satisfied per structure (isolated/in situ)']=validation.pop('Percent of excluded volume restraints satisfied per structure')


print('---------------')
print(samples)

samples['Number of models after equilibration']= '300000/300000'
samples['Number of models that satisfy the input information']= '92281/51183'
samples['Number of structures in samples A/B']= '52831,39450/18722,32461' 
samples['p-value of non-parametric Kolmogorov-Smirnov two-sample test']= '0.446/0.34'
samples['Kolmogorov-Smirnov two-sample test statistic, D']= '0.0,1E-16'

samples['Number of models after equilibration (isolated/in situ)']=samples.pop('Number of models after equilibration')
samples['Number of models that satisfy the input information (isolated/in situ)']=samples.pop('Number of models that satisfy the input information')
samples['Number of structures in samples A,B (isolated/in situ)']=samples.pop('Number of structures in samples A/B')
samples['p-value of non-parametric Kolmogorov-Smirnov two-sample test (isolated/in situ)']=samples.pop('p-value of non-parametric Kolmogorov-Smirnov two-sample test')
samples['Kolmogorov-Smirnov two-sample test statistic (D) (isolated/in situ)']=samples.pop('Kolmogorov-Smirnov two-sample test statistic, D')


###################
print('-------------')
print(clustering)


clustering['Sampling precision']='14.66/12.6 \\AA'
clustering["Homogeneity of proportions $\\chi^2$ test (p-value)/Cramer???s V value"]="1.000 (0.000)/ 1.000 (0.000) (thresholds: p-value$>$0.05 OR Cramer's V$<$0.1)"
#"1.000 (0.000)/ 1.000 (0.000)(thresholds: p-value$>$0.05 OR Cramer's V$<$0.1)"]
clustering['Number of clusters']='1/1'
clustering['Cluster populations']='Cluster 1: 99/98 \%'
clustering['Cluster precisions']='Cluster 1: 15.42/10.1 \\AA'
clustering['Average cross-correlation between localization probability densities of samples A and B']='Cluster 1: 0.84/0.92'

clustering['Sampling precision (isolated/in situ)']=clustering.pop('Sampling precision')
clustering['Homogeneity of proportions $\\chi^2$ test p-value (Cramer???s V value) (isolated/in situ)']=clustering.pop('Homogeneity of proportions $\\chi^2$ test (p-value)/Cramer???s V value')
clustering['Number of clusters (isolated/in situ)']=clustering.pop('Number of clusters')
clustering['Cluster populations (isolated/in situ)']=clustering.pop('Cluster populations')
clustering['Cluster precisions (isolated/in situ)']=clustering.pop('Cluster precisions')
clustering['Average cross-correlation between localization probability densities of samples A and B (isolated/in situ)']=clustering.pop('Average cross-correlation between localization probability densities of samples A and B')


software['Modeling scripts'] = ['https://github.com/integrativemodeling/NPC\_v3.0']
software['Homology detection and structure prediction'] = ['HHPred, version 2.0.16']
software['Visualization and plotting'] = ['UCSF Chimera', 'Matplotlib, version 3.0.3 ']

################################################
# Convert ordered dictionaries 
# into lists
################################################
input_information_list = dict_to_list(input_information)
representation_list = dict_to_list(representation)
sampling_list = dict_to_list(sampling)
enumeration_list = dict_to_list(enumeration)
validation_threading_list = dict_to_list(validation_threading)
samples_list = dict_to_list(samples)
clustering_list = dict_to_list(clustering)
validation_list = dict_to_list(validation)
software_list = dict_to_list(software)
data_availability_list = dict_to_list(data_availability)

print(enumeration_list)


print(sampling_list)


################################################
# Compile all information
# 
################################################
variable_dict = {'complex': 'NPC',
                 'number':1,
                 'input_information': input_information_list, 
                 'representation': representation_list,
                 'enumeration': enumeration_list,
                 'sampling': sampling_list,
                 'validation_threading':validation_threading_list,
                 'samples': samples_list,
                 'clustering':clustering_list,
                 'validation':validation_list,
                 #'benchmark':benchmark_list,
                 'software':software_list,
                 'data':data_availability_list}

################################################
# Generate tex, pdf file
################################################
template = utils.get_template('../utils/SI_template.tex')
utils.compile_pdf_from_template(template, variable_dict, './table_SI_NPC_v3.pdf')

exit()
