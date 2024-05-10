#!/bin/bash
#$ -cwd

rm script_SLIM3_2R_total.sh.*

if [ $# -ne 3 ]
then
	echo "Usage: $0 <slim3INPUTfile> <REPS> <n>"
	exit 1
fi

#Set arguments
INPUT=$1
REPS=$2
n=$3

#Working directory
WDIR=$PWD

#Output directory
if [ -d "$WDIR/OUTPUT_SLIM_$n" ]
then
rm -r $WDIR/OUTPUT_SLIM_$n
fi

mkdir -p $WDIR/OUTPUT_SLIM_$n
mkdir -p /state/partition1/SLIM$n/$SLURM_JOBID/

###################### TRANSFER TO state/partition1 #########################
cp SNP_BP_SLIM3_2 /state/partition1/SLIM$n/$SLURM_JOBID/
cp SIMOVERDOM /state/partition1/SLIM$n/$SLURM_JOBID/
cp AVERAGE /state/partition1/SLIM$n/$SLURM_JOBID/
cp $INPUT /state/partition1/SLIM$n/$SLURM_JOBID/
cp recombination-file-* /state/partition1/SLIM$n/$SLURM_JOBID/
cp genes* /state/partition1/SLIM$n/$SLURM_JOBID/
cp B-* /state/partition1/SLIM$n/$SLURM_JOBID/
cp shell_calculations /state/partition1/SLIM$n/$SLURM_JOBID/
cp plink /state/partition1/SLIM$n/$SLURM_JOBID/
cp vcftools /state/partition1/SLIM$n/$SLURM_JOBID/

touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`
cd /state/partition1/SLIM$n/$SLURM_JOBID

################################ REPLICATES #############################
for ((r=1; r<=$REPS; r++))
do

###################### run SLIM3 #########################
module load SLiM/3.3.2

START=$(date +%s)
slim $INPUT > slimout$r
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Slim3 took 		$DIFF seconds" >> timefile

cp collect_par_mutations.txt collect_par_mutations$r.txt

###################### SNP_BP_SLIM3_2 #########################
num=$RANDOM
echo "$num" > seedfile

cp slimout$r slimout

START=$(date +%s)
./SNP_BP_SLIM3_2<<@
-99
1000	N
@
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "SNP_BP took 	$DIFF seconds" >> timefile

cp dataBP.ped dataBP$r.ped
cp dataBP.map dataBP$r.map
cp list_allsnps list_allsnps$r
cp list_qtls list_qtls$r

###################### SIMOVERDOM ########################
num=$RANDOM
echo "$num" > seedfile

START=$(date +%s)
./SIMOVERDOM<<@
1
-99
1000    NINDNP (max 10000)
50     NIND (max 10000)
100             classes
21100000          nwindows
100000          window_size
100              Replicates
@
DIFF=$(( $END - $START ))
echo "SIMOVERDOM took 		$DIFF seconds" >> timefile

cp genfile.dat genfile$r.dat
cp table.txt table$r.txt
rm table.txt

######################## AVERAGE #########################
START=$(date +%s)
./AVERAGE<<@
$r    NREPS (max 100)
211   NWINDOWS
@
DIFF=$(( $END - $START ))
echo "AVERAGE took 		$DIFF seconds" >> timefile

cp timefile timefile$r

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r /state/partition1/SLIM$n/$SLURM_JOBID/dataBP$r.ped $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/dataBP$r.map $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/list_allsnps$r $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/list_qtls$r $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/slimout$r $WDIR/OUTPUT_SLIM_$n
#cp -r /state/partition1/SLIM$n/$SLURM_JOBID/collect_par_mutations$r.txt $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/genfile$r.dat $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/table$r.txt $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/tableavg.txt $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/se.txt $WDIR/OUTPUT_SLIM_$n
#cp -r /state/partition1/SLIM$n/$SLURM_JOBID/data.map $WDIR/OUTPUT_SLIM_$n
#cp -r /state/partition1/SLIM$n/$SLURM_JOBID/data.ped $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/timefile$r $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/tablepi* $WDIR/OUTPUT_SLIM_$n
cp -r /state/partition1/SLIM$n/$SLURM_JOBID/tableDT* $WDIR/OUTPUT_SLIM_$n

rm collect_par_mutations*
rm slimout*
rm data*
rm list*

done

######################## state/partition1 CLEANING #########################

rm -r /state/partition1/SLIM$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*