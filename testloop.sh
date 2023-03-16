#!/bin/bash
for i in {0..8}
do
#   while [[ $numproc -le 8]]
    for numproc in 1 2 4
    do
        echo "testcase : $i"
        echo "num procs : $numproc"
        time mpirun -np $numproc ./a2 --taskid=1 --inputpath=./../testcases/A2/test$i/test-input-$i.gra --headerpath=./../testcases/A2/test$i/test-header-$i.dat --outputpath=outputfile.txt --verbose=0 --startk=1 --endk=30 --p=2
    #   numproc = $((numproc * 2))
    done
done
