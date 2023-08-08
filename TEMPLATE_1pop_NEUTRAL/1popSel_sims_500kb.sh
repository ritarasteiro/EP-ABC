#!/bin/bash
#PBS -N CHANGEJOBNAME
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -j oe
#PBS -o out/
#PBS -l walltime=CHANGE_WALLTIME
# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m a
#PBS -M ar17162@bristol.ac.uk
#PBS -J 1-CHANGE_nbARRAYS

cd $PBS_O_WORKDIR

##module for R
module load lang/r/3.6.1
##modules for msprime and slim
module add tools/cmake/3.14.2
module add lang/gcc/9.1.0
module load lang/python/anaconda/3.7-2019.03


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CHANGES IN THIS SCRIPT 8 OCT2019
# Copy infiles to $TMPDIR : to decrease computation time
#CHANGES IN THIS SCRIPT 9 OCT2019
# removed creation of folder named after array and instead
# changed name of output files to have $PBS_ARRAY_INDEX: avoids overwriting of different outputs created in different arrays

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# record some potentially useful details about the job:
echo "The Array ID is: $PBS_ARRAY_INDEX"
echo Running on host `hostname`
echo Start time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This jobs runs on the following machines:

nbITER=CHANGE_NITER
iter=$(cat iter_nb) # keep track of the iteration number

if [ $iter -gt $nbITER ]
then
	exit

else

	cd iter$iter
	
	# nb of sims per array
	nSIMS=CHANGE_MAX_ARRAY
	ITNB=$(awk '{print $1}' iter_nb)

	#calculate sample size (haploid) # CHANGE value depending of the species
	SAMPLE_IND=CHANGE_IND
	SAMPLE=$(r -e "cat(as.integer($SAMPLE_IND*2),sep='\n')")


	

#	copy infiles to TMPDIR
	cp -r ../files_to_run/. $TMPDIR/
	cp  params.sims$ITNB.txt $TMPDIR/

#	# set folder for simulation
	FOLDER=$PBS_ARRAY_INDEX
	mkdir $FOLDER


	#############################################################################
	# Establish start and end index of sims per array
	START=$(r -e "cat(as.integer(seq(1,CHANGE_TOTALNBSIMS,CHANGE_MAX_ARRAY)[$PBS_ARRAY_INDEX]),sep='\n')")
	# end sims
	END=$(r -e "cat(as.integer($START+$nSIMS-1),sep='\n')")

	for (( sim=$START; sim<=$END; sim++ ))
	do

		echo SIMULATION  $sim
		cd $FOLDER

		##################################################################################
		
		sed -n $sim"p" $TMPDIR/params.sims$ITNB.txt > temp$PBS_ARRAY_INDEX

		#Ne 
		N0=$(awk '{print $1}' temp$PBS_ARRAY_INDEX)
		N1=$(awk '{print $2}' temp$PBS_ARRAY_INDEX)
		N2=$(awk '{print $3}' temp$PBS_ARRAY_INDEX)
		N3=$(awk '{print $4}' temp$PBS_ARRAY_INDEX)
		N4=$(awk '{print $5}' temp$PBS_ARRAY_INDEX)
		N5=$(awk '{print $6}' temp$PBS_ARRAY_INDEX)
		N6=$(awk '{print $7}' temp$PBS_ARRAY_INDEX)
		N7=$(awk '{print $8}' temp$PBS_ARRAY_INDEX)
		N8=$(awk '{print $9}' temp$PBS_ARRAY_INDEX)
		N9=$(awk '{print $10}' temp$PBS_ARRAY_INDEX)
		N10=$(awk '{print $11}' temp$PBS_ARRAY_INDEX)
		Nanc=$(awk '{print $12}' temp$PBS_ARRAY_INDEX)

		#time parameters
		T1=$(awk '{print $13}' temp$PBS_ARRAY_INDEX)
		T2=$(awk '{print $14}' temp$PBS_ARRAY_INDEX)
		T3=$(awk '{print $15}' temp$PBS_ARRAY_INDEX)
		T4=$(awk '{print $16}' temp$PBS_ARRAY_INDEX)
		T5=$(awk '{print $17}' temp$PBS_ARRAY_INDEX) 
		T6=$(awk '{print $18}' temp$PBS_ARRAY_INDEX)
		T7=$(awk '{print $19}' temp$PBS_ARRAY_INDEX)
		T8=$(awk '{print $20}' temp$PBS_ARRAY_INDEX) 
		T9=$(awk '{print $21}' temp$PBS_ARRAY_INDEX)
		T10=$(awk '{print $22}' temp$PBS_ARRAY_INDEX)
		T11=$(awk '{print $23}' temp$PBS_ARRAY_INDEX) 
#		#mutation and recombination rates
		MUT=$(awk '{print $24}' temp$PBS_ARRAY_INDEX) 
		REC=$(awk '{print $25}' temp$PBS_ARRAY_INDEX) 

		rm temp$PBS_ARRAY_INDEX

		#number of bases
		NBASES="50e+04" #500kb windows

		# to avoid N smaller than sample size OR very large N to be simulated in slim (very slow) 
		if [ $N0 -gt 500000 -o  $N1 -gt 500000 -o  $N2 -gt 500000 -o  $N3 -gt 500000 -o $N4 -gt 500000 -o $N5 -gt 500000 -o $N6 -gt 500000 -o  $N7 -gt 500000 -o  $N8 -gt 500000 -o $N9 -gt 500000 -o $N10 -gt 500000 -o $Nanc -gt 500000 ]
		then
		echo 'NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA' >> ../summary_stats_$PBS_ARRAY_INDEX.txt
		echo $sim >> ../param_index_$PBS_ARRAY_INDEX # keep the indices of the combination of parameters simulated (check for when nodes are not working properly and freeze)
		echo "one N is too large"		
		cd ../
	
		else

			##############################
			#change msprime infile
			##############################
			#open .py file and replace parameters
			sed -s "s/\<NBASES\>/${NBASES}/g;s/MUT/${MUT}/g;s/REC/${REC}/g;s/\<N0\>/${N0}/g;s/\<N1\>/${N1}/g;s/\<N2\>/${N2}/g;s/\<N3\>/${N3}/g;s/\<N4\>/${N4}/g;s/\<N5\>/${N5}/g;s/\<N6\>/${N6}/g;s/\<N7\>/${N7}/g;s/\<N8\>/${N8}/g;s/\<N9\>/${N9}/g;s/\<N10\>/${N10}/g;s/\<Nanc\>/${Nanc}/g;s/\<T1\>/${T1}/g;s/\<T2\>/${T2}/g;s/\<T3\>/${T3}/g;s/\<T4\>/${T4}/g;s/\<T5\>/${T5}/g;s/\<T6\>/${T6}/g;s/\<T7\>/${T7}/g;s/\<T8\>/${T8}/g;s/\<T9\>/${T9}/g;s/\<T10\>/${T10}/g;s/\<T11\>/${T11}/g;s/\<SAMPLE\>/${SAMPLE}/g" $TMPDIR/model0.py >  GrM.infile.py
			
			##############################
			#RUN MSPRIME
			##############################
	 		python GrM.infile.py
			
			# if position_slim is empty
			if [[ -z $(grep '[^[:space:]]' position) ]]
			then
			echo '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA' >> ../summary_stats_$PBS_ARRAY_INDEX.txt
			echo $sim >> ../param_index_$PBS_ARRAY_INDEX # keep the indices of the combination of parameters simulated (check for when nodes are not working properly and freeze)
			echo "msprime output monomorphic"
			cd ../
			
			else
				#CREATE a MS_type input file
				s=$(wc -l < position)
				position=$(r -e "options(scipen=999);cat(scan(argv[1],sep=' ')/$NBASES)" position )
				echo "//" > ms_out
				echo "segsites:" $s >> ms_out
				echo "positions:" $position >> ms_out
				paste loci>> ms_out
				rm loci position
				cd ../
				Rscript ../files_to_run/run_sumstatGrM.r $PBS_ARRAY_INDEX/ms_out $SAMPLE $NBASES $PBS_ARRAY_INDEX
				
				#check if R crashes (bug within Rcpp)
				Rstatus=$(echo $?) #status of previous shell command (Rscript): 0 (exit without errors) or 1 (exit with errors) 
				if [ $Rstatus -eq 0 ]; then
 				echo $sim >> param_index_$PBS_ARRAY_INDEX # keep the indices of the combination of parameters simulated (check for when nodes are not working properly and freeze)

				#if R crashes
				else
				echo 'NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA' >> ../summary_stats_$PBS_ARRAY_INDEX.txt
				echo $sim >> param_index_$PBS_ARRAY_INDEX # keep the indices of the combination of parameters simulated (check for when nodes are not working properly and freeze)
				echo  "RCPP crashed"
				fi
					
			fi
		fi
	
	
	done
	rm -r $PBS_ARRAY_INDEX
	cd ../

fi
cd $PBS_O_WORKDIR

echo End time is `date`
