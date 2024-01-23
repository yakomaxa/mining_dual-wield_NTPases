#!/bin/sh
source /data/group1/z44714a/gromacs/bin/bin/GMXRC.bash
module unload cuda/11.2.1
module load cuda/11.6.2
module unload gcc
module load gcc/10.3.0

# generaete tpr for nvt-eq
gmx grompp -f ./mdps/nvt.mdp -c short_nvted_${1}.gro  -t short_nvted_${1}.cpt -p topol_${1}.top -o nvt_${1}.tpr -n index.ndx -r short_nvted_${1}.gro
# run nvt-eq
gmx mdrun -v -deffnm nvted_${1} -s nvt_${1}.tpr -nb gpu

# generaete tpr for npt-eq
gmx grompp -f ./mdps/npt.mdp -c nvted_${1}.gro -p topol_${1}.top -o npt_${1}.tpr -t nvted_${1}.cpt -n index.ndx -r nvted_${1}.gro
# run npt-eq
gmx mdrun -v -deffnm npted_${1} -s npt_${1}.tpr -nb gpu

# generaete tpr for production
gmx grompp -f ./mdps/short_md.mdp -c npted_${1}.gro -t npted_${1}.cpt -p topol_${1}.top -o short_produce_${1}.tpr -n index.ndx 
# production run!
gmx mdrun  -deffnm short_produced_${1}  -s short_produce_${1}.tpr -nb gpu
