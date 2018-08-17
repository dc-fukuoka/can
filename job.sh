#!/bin/sh
#$ -cwd
#$ -l f_node=64
#$ -l h_rt=1:0:0
#$ -p -3

. /etc/profile.d/modules.sh
module purge
module load cuda pgi openmpi
module load intel
module list 2>&1

nnodes=64

mpirun -npernode 4 -np $((4*nnodes)) ./can_acc
