#!/bin/csh

#$ -M anematba@nd.edu	 # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  gpu 	 # Specify queue
#$ -l gpu_card=1
#$ -pe smp 4 
#$ -N  run_May11	 # Specify job name


module load gcc/
module load cuda
echo -n "It is currently: ";date
echo -n "I am logged on as ";who am i
echo -n "This computer is called ";hostname
echo -n "I am currently in the directory ";pwd
#setenv PATH /afs/crc.nd.edu/user/a/anematba/Public/2015/Oct/11th/SceCells/bin:$PATH
./bin/runDiscSimulation_M -slurm N01_2
