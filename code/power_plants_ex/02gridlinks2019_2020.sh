#!/bin/sh

## Give your job a name to distinguish it from other jobs you run.
#SBATCH --job-name=02gridlinks19-20

## General partitions: all-HiPri, bigmem-HiPri   --   (12 hour limit)
##                     all-LoPri, bigmem-LoPri, gpuq  (5 days limit)
## Restricted: CDS_q, CS_q, STATS_q, HH_q, GA_q, ES_q, COS_q  (10 day limit)
#SBATCH --partition=normal #all-LoPri

## Separate output and error messages into 2 files.
## NOTE: %u=userID, %x=jobName, %N=nodeID, %j=jobID, %A=arrayID, %a=arrayTaskID
#SBATCH --output=/scratch/%u/logs/%x-%a.out  # Output file
#SBATCH --error=/scratch/%u/logs/%x-%a.err   # Error file

## Slurm can send you updates via email
#SBATCH --mail-type=BEGIN,END,FAIL         # ALL,NONE,BEGIN,END,FAIL,REQUEUE,..
#SBATCH --mail-user=lhennem@gmu.edu     # Put your GMU email address here

## Specify how much memory your job needs. (2G is the default)
#SBATCH --mem-per-cpu=10G        # Total memory needed per task (units: K,M,G,T)

## Specify how much time your job needs. (default: see partition above)
#SBATCH --time=01-06:00  # Total time needed for job: Days-Hours:Minutes

## create an array of jobs
#SBATCH -a 1-4

## select a specific node
#SBATCH --nodelist=hop070

## Load the relevant modules needed for the job
module load r-disperseR/0.1.0

## Run your program or script
Rscript /home/lhennem/runscripts/rundisperseR/02gridlinks2019_2020.R 

#source( "/n/home03/lhenneman/runscripts/rundisperseR/01disperseR2016_2018.R")

