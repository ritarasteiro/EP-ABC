#!/bin/bash
#PBS -N mergeABC
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -j oe
#PBS -o out/
#PBS -l walltime=CHANGE_WALLTIME2
# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m a
#PBS -M ar17162@bristol.ac.uk


cd $PBS_O_WORKDIR

##module for R
module load lang/r/3.6.1


# record some potentially useful details about the job:
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This jobs runs on the following machines:

#arguments for epABC
MODEL=CHANGEMODEL # CHANGE pop name
GEN=CHANGE_GEN #nb of generations per year # CHANGE value depending of the species
nbITER=CHANGE_NITER
nbSIMS=50000
TOL=0.05
ESS_nbSIMS_RATIO=0.3
SCALE=CHANGE_SCALE
iter=$(cat iter_nb) # keep track of the iteration number
NLOCI=CHANGE_LOCI
echo $iter


if [ ! -f iter_nb ] 
then
#merge projection sims
cd PROJECTION
find ./ -name "summary_stats_*.txt" | sort -V | xargs cat >> SS_proj.txt
find ./ -name "par.sims.proj*" | sort -V | xargs cat >> params.sims.proj.txt
find ./ -name "par.proj*.txt" | sort -V | xargs cat >> params.proj.txt

nbSIMS_proj=100000

wc SS_proj.txt
wc params.proj.txt
	s=$(wc -l < SS_proj.txt)
	p=$(wc -l < params.proj.txt)

	if [ $s -eq $nbSIMS_proj -a  $p -eq $nbSIMS_proj ]
	then
	rm summary_stats_*.txt
	rm param_index_*
	rm par.proj*.txt
	rm par.sims.proj*.txt
	cd ..
	# create infiles for iteration 1 
	Rscript iter1_params.R $nbSIMS $SCALE $GEN $NLOCI
	else
		exit
	fi

#end of iterations
elif [ $iter -gt $nbITER ]
then
exit

#iterations>=1
else
cd iter$iter
#pwd
	#MERGE files in numerical order
	find ./ -name "summary_stats_*.txt" | sort -V | xargs cat >> SS$iter.txt
	find ./ -name "param_index_*" | sort -V | xargs cat >> params_indices$iter.txt
	wc SS$iter.txt
	wc params_indices$iter.txt
	s=$(wc -l < SS$iter.txt)
	p=$(wc -l < params_indices$iter.txt)
	
	if [ $s -eq $nbSIMS -a  $p -eq $nbSIMS ]
	then
	rm summary_stats_*.txt
	rm param_index_*
	
	#RUN EP-ABC: runs statistical analyses and creates new inputs for the simulator. 
	#Our method recicles simulations through important sampling, thus we do not need to simulate at each iteration 
	cd ../
#pwd
	Rscript epABC_proj.R $MODEL $nbITER $nbSIMS $TOL $ESS_nbSIMS_RATIO $SCALE $GEN 
	
	else
		#Checking for errors in runs (frozen nodes, etc)
		echo 'iteration ' $iter >>../iter$iter.$PBS_JOBID
		for i in {1..CHANGE_nbARRAYS}
		do
		ss=$(wc -l < summary_stats_$i.txt)
		pa=$(wc -l < param_index_$i)
		nbSIMS_array=CHANGE_MAX_ARRAY
		if [ $ss -eq $nbSIMS_array -a  $pa -eq $nbSIMS_array ]
		then
		echo ARRAY $i OK >> ../iter$iter.$PBS_JOBID
		else
		echo ARRAY $i ERROR SS = $ss param = $pa >>../iter$iter.$PBS_JOBID
		fi
		done

	rm summary_stats_*.txt
	rm param_index_*
	rm SS$iter.txt
	rm params_indices$iter.txt

	fi


fi

cd $PBS_O_WORKDIR
echo End time is `date`
