#!/bin/sh
# header to run at Golgi server
# #PBS -V
# #PBS -l nodes=1:ppn=12
# #cd $PBS_O_WORKDIR
#source /usr/local/gromacs/bin/GMXRC

# source Gromacs setting file for elm0
#source /local/apl/lx/gromacs2020.6-CUDA/bin/GMXRC.zsh
source /data/group1/z44714a/gromacs/bin/bin/GMXRC.bash
source /data/group1/z44714a/gromacs/bin/bin/GMXRC.bash
# make topology using force field in the directory (amberff99 )
gmx pdb2gmx -f ${1} -o processed_${1}.gro -water tip3p -p topol_${1}.top -ignh -ter < ./dummy/1-1-1.txt

# make box
gmx editconf -f processed_${1}.gro -o newbox_${1}.gro -c -d 1.3 -bt dodecahedron
#gmx editconf -f processed_${1}.gro -o newbox_${1}.gro -c -bt cubic -box 6

# generate tpr (input for mdrun) for minimzation in vacuo
gmx grompp -f ./mdps/minim.mdp -c newbox_${1}.gro -p topol_${1}.top -o vacuoem_${1}.tpr -maxwarn 10

# run in-vacuo minimization
gmx mdrun -deffnm vacuoemed_${1} -nb gpu -s vacuoem_${1}.tpr

# solvate (place water model)
gmx solvate -cp vacuoemed_${1}.gro -cs spc216.gro -o solv_${1}.gro -p topol_${1}.top 

# make topology for ions
gmx grompp -f ./mdps/ions.mdp -c solv_${1}.gro -p topol_${1}.top -o ions_${1}.tpr -maxwarn 10

# replace some water by ions
gmx genion -s ions_${1}.tpr -o solv_ions_${1}.gro -p topol_${1}.top -pname NA -nname CL  -neutral  -conc 0.1 -seed 1 < ./dummy/18.txt

# generate tpr for  minimzation with water and ions
gmx grompp -f ./mdps/minim.mdp -c solv_ions_${1}.gro -p topol_${1}.top -o em_${1}.tpr 

# perform minimzation with water and ions
gmx mdrun  -deffnm emed_${1} -nb gpu -s em_${1}.tpr 
