#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=6:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=50G   # memory per CPU core
#SBATCH -J "all"   # job name
##SBATCH --qos=test

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

#python fasterCombineFiles.py bacteria_files_longestIsoform/
#time python cUsageNew.py -id all_filtered/ -o all_out_mine2
#time python cUsageNew.py -id bacteria_filtered/ -o bacteria_distMatrix
#time python fastNJ.py bacteria_distMatrix_revised bacteria_trees.nwk
#time python extractTNTByGeneNames.py -id bacteria_filtered/ -o geneOnly_bacteria/ -g
#time python extractTNTByGeneNames.py -id all_filtered/ -o geneOnly_all/ -g
#cp -r geneOnly_* ../forLauren/

#time python cam.py -aa -id ../../all_filtered_noPartial/ >all_aa_matrix
echo "Bacteria Neighbor joining"
#time python ../probablyDontNeed/fastNJ.py aa_all_matrix aa_all_tree
time python ../probablyDontNeed/fastNJ.py aa_bacteria_matrix aa_bacteria_tree2
#time python cam.py -aa -id ../../bacteria_filtered >aa_all_matrix
#time python cam.py  -id ../../all_filtered >ha2



#time python makeNewickFromDistBioPythonNJ.py all_out_mine_ordered >all_nj.nwk
#time python makeNewickFromDistBioPython.py all_out_mine_ordered >all_upgma.nwk
#time python makeNewickFromDistBioPython.py bacteria_distMatrix >bacteria_upgma.nwk
#time python makeNewickFromDistBioPythonNJ.py bacteria_distMatrix >bacteria_nj.nwk
#time python makeNewickFromDistBioPython.py viral_distMatrix >viral_upgma.nwk
#time python makeNewickFromDistBioPythonNJ.py viral_distMatrix >viral_nj.nwk




#time python nj.py all_out_mine_ordered >all_out_mine_nj

#ls bacteria_files/ | xargs -P 16 -n 1 -I {} python findLongestIsoform.py bacteria_files/'{}' bacteria_files_longestIsoform/'{}'

#tar -xf bacteria.tar
#python discretecounts.py $1
#python codonfrequency.py gff3_parser/filteredGenes_noException/$1 gff3_parser/codon_bias_noException/$1
#chmod -R 777 ALL_genomes
#python calculatecub.py $1
#python comparevirusprot.py $1
#./doIt.sh $1

#input="CDS/"$1
#output="ONLY_CDS/"$1$1
#out2="ONLY_CDS/"$1


#awk '!a[$0]++' $input >$output
#python eraseHeaders.py $output $out2
#rm $output
#python gff3_parser.py ../../ALL_genomes/H_sapiens/GFF/ref_GRCh38_scaffolds.gff3.gz H_sapiens
#python gff3_parser.py $1 $2

#./doIt.sh $1


exit 0
