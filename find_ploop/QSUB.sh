# fullpath is a fullpath for splited index file for foldcomp. 
# node is the node name to run the chunk
# you need this kind of file shown below as input
# /your_path_to_chunked_index/xaa node01
# /your_path_to_chunked_index/xab node02
# /your_path_to_chunked_index/xac node03
# /your_path_to_chunked_index/xad node04

while read fullpath node
do 
    name=`basename ${fullpath}`
    cp -r template dir_${name} 
    cd dir_${name} 
    qsub -q all.q@${node} -N fragQ_${name} run-frag.sh ${fullpath}
    cd ../ 
done < ${1}
