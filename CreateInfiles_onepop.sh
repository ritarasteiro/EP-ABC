#bash
#create replicate folders
# change the parameters for the runs
MODEL=Model1
JOBNAME=Model1
TEMPLATE=TEMPLATE_1pop

NLOCI=100

LENGTH=500
SCALE=2

SAMPLE_IND=30
GEN=3
nbSIMS=50000
nbARRAYS_JOB=200

NBASES=$LENGTH"e+03" #in KB
SS=46  #Total nb of ss
NA=$(r -e "cat(as.integer($SS-20))") #nb of ss with na when pop is monomorphic (this corresponds to r2 and roh. All the other 20 ss will be equal to 1)

nbITER=1000
WALLTIME1="16:00:00"
WALLTIME2="3:00:00"
MAX_ARRAY=$(r -e "cat(as.integer($nbSIMS/$nbARRAYS_JOB),sep='\n')") #nb sims per array


# Do not change code below
#with parallelABC
mkdir $MODEL
cd $MODEL
#cp -r ../PROJECTION$LENGTH$EVOL/PROJECTION/ ./
for i in {1..5}
do
mkdir run$i
mkdir run$i/out
cp -r ../$TEMPLATE/. run$i/
sed -s "s/\<CHANGEJOBNAME\>/$JOBNAME.$i/g;s/\<CHANGE_IND\>/$SAMPLE_IND/g;s/\<CHANGE_NITER\>/$nbITER/g;s/\<CHANGE_WALLTIME\>/$WALLTIME1/g;s/\<CHANGE_MAX_ARRAY\>/$MAX_ARRAY/g;s/\<CHANGE_TOTALNBSIMS\>/$nbSIMS/g;s/\<CHANGE_nbARRAYS\>/$nbARRAYS_JOB/g;s/\<CHANGE_NBASES\>/$NBASES/g;s/\<CHANGE_NBSS\>/$SS/g;s/\<CHANGE_NA\>/$NA/g"   run$i/1popSel_sims.sh >temp
mv temp  run$i/1popSel_sims.sh
sed -s "s/\<CHANGE_GEN\>/$GEN/g;s/\<CHANGE_LOCI\>/$NLOCI/g;s/\<CHANGE_SCALE\>/$SCALE/g;s/\<CHANGEMODEL\>/$MODEL/g;s/\<CHANGE_NITER\>/$nbITER/g;s/\<CHANGE_WALLTIME2\>/$WALLTIME2/g;s/\<CHANGE_MAX_ARRAY\>/$MAX_ARRAY/g;s/\<CHANGE_nbARRAYS\>/$nbARRAYS_JOB/g"  run$i/merge_EPabc.sh > temp
mv temp  run$i/merge_EPabc.sh


cp ../$MODEL'_'targetSS.txt run$i/

cd run$i
Rscript iter1_params.R $nbSIMS $SCALE $GEN $NLOCI $MODEL
cd ..
rm temp
done 
cd ..




