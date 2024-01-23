#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -pe make 4

## -pe make 4 is declartion of number of core used 

# export shared library for gromacs configured for AMD CPU
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/sakuma/.gromacs-2016.3_gcc5.3_AMD/libcopy/usr/local/lib64/:/home/sakuma/.gromacs-2016.3_gcc5.3_AMD/libcopy/usr/lib64/

# source gromacs setting files for AMD CPU ( elm1 elm2 elm3 elm4 )
#source /home/sakuma/.gromacs-2016.3_gcc5.3_AMD/bin/GMXRC.zsh
source /data/group1/z44714a/gromacs/bin/bin/GMXRC.bash

# generaete tpr for nvt-eq
gmx grompp -f ./mdps/nvt.mdp -c short_nvted_${1}.gro  -t short_nvted_${1}.cpt -p topol_${1}.top -o nvt_${1}.tpr -n index.ndx
# run nvt-eq
gmx mdrun -v -deffnm nvted_${1} -s nvt_${1}.tpr -ntomp 8  -ntmpi 8

# generaete tpr for npt-eq
gmx grompp -f ./mdps/npt.mdp -c nvted_${1}.gro -p topol_${1}.top -o npt_${1}.tpr -t nvted_${1}.cpt -n index.ndx
# run npt-eq
gmx mdrun -v -deffnm npted_${1} -s npt_${1}.tpr -ntomp 8  -ntmpi 8

# generaete tpr for production
gmx grompp -f ./mdps/short_md.mdp -c npted_${1}.gro -t npted_${1}.cpt -p topol_${1}.top -o short_produce_${1}.tpr -n index.ndx
# production run!
gmx mdrun  -deffnm short_produced_${1}  -s short_produce_${1}.tpr -ntomp 8  -ntmpi 8 -deffnm
