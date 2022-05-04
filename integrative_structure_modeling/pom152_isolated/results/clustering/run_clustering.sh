module load python3/pyrmsd
module load python3/numpy
module load imp

export top_dir=/home/ignacia/Research/yeast_NPC/membrane_ring/manual_fitting_manual_noshuffle_v3/
export analys_dir=$top_dir/analys
#export mod_dir=$analys_dir/GSMs_cl0/
export name=pom152
export cluster=0

#cp $top_dir/density.txt $mod_dir
#cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
#cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

cp $top_dir/density.txt .
#cp $mod_dir/sample_A.txt Scores_A.txt
#cp $mod_dir/sample_B.txt Scores_B.txt

#ls -lta $analys_dir/GSMs_cl${cluster}/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster${cluster}_detailed.dat
#ls -lta $analys_dir/GSMs_cl${cluster}/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster${cluster}_detailed.dat

nohup python3 ~/SOFTW/imp-sampcon-0320/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $analys_dir --mode cuda --cores 4 --align --density density.txt --gridsize 3.0 --gnuplot --scoreA $analys_dir/A_models_clust${cluster}.txt --scoreB $analys_dir/B_models_clust${cluster}.txt  --rmfA $analys_dir/A_models_clust${cluster}.rmf3 --rmfB $analys_dir/B_models_clust${cluster}.rmf3 -dt 30. > clustering.log &
