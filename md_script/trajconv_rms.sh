gmx trjconv -s short_produce_${1}.pdb.tpr  -f npted_${1}.pdb.gro  -o ${1}_noPBC.gro -pbc mol -center < ../1.txt
gmx trjconv -s short_produce_${1}.pdb.tpr  -f short_produced_${1}.pdb.xtc  -o ${1}_noPBC.xtc -pbc mol -center < ../1.txt
gmx rms -s em_${1}.pdb.tpr -f ${1}_noPBC.xtc -o rmsd_xtal_${2}.xvg -tu ns < ../3.txt
