source /data/group1/z44714a/gromacs/bin/bin/GMXRC.bash
module unload cuda/11.2.1
module load cuda/11.6.2
module unload gcc
module load gcc/10.3.0
gmx trjconv -s short_produce_${1}.tpr -f short_produced_${1}.xtc -o trajectory_wo_pbc.xtc -pbc mol -ur compact -center < ./1-1.txt
 
