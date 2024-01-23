for i in {01..01} 
do 
 mkdir S${i}
 cd S${i} 
 for j in {01..20} 
  do cp -r ../template traj${j} 
  cp ../fullatom/${i}.pdb traj${j} 
   cd traj${j} 
   bash ./topol_box_vacuo_solv_ion_minim.sh ${i}.pdb 
   bash ./short_nvt.sh ${i}.pdb 
   cat ./submit.sh | sed "s/INPUTST/${i}.pdb/g" > tmp ; mv tmp ./submit.sh 

   cd ../ 
done 
cd ../ 
done
