#/bin/bash

echo '' > results.txt;

while read -r line;
do
    echo $line;
    echo $line > matfileReal;
    echo "##" >> matfileReal;
    #~ cp ./inputs_vbarms ./inputs;
    #~ mpiexec -n 3 ./Bdd-mtx-dse.ex >> results.txt;
    cp ./inputs_arms ./inputs;
    #~ ./dd-mtx-dse.ex >> results.txt;
    mpiexec -n 3 ./dd-mtx-dse.ex >> results.txt;
done < matrices.txt
