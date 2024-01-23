#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V

# Run quite short nvt simulation (with restraint) to check that model do not explode. If the model explodes, it is not worth running longer md.

# source setting file for elm0
#source /usr/local/gromacs-5.1.2_gcc5.3/bin/GMXRC.zsh
#source /local/apl/lx/gromacs2020.6-CUDA/bin/GMXRC.zsh
source /data/group1/z44714a/gromacs/bin/bin/GMXRC.bash
# generate index specifying protein atoms
gmx make_ndx -o index.ndx -f emed_${1}.gro < ./dummy/1.txt

# generate tpr for short nvt run
gmx grompp -f ./mdps/short_nvt.mdp -c emed_${1}.gro -p topol_${1}.top -o short_nvt_${1}.tpr  -n index.ndx -r emed_${1}.gro

# run test md
gmx mdrun -deffnm short_nvted_${1} -nb gpu -s short_nvt_${1}.tpr 

# If you do not find the model explodes, go next for longer nvt, npt equilibration and production run !
