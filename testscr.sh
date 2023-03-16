time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test0/test-input-0.gra --headerpath=./../testcases/A2/test0/test-header-0.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=10 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test0/output0_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test1/test-input-1.gra --headerpath=./../testcases/A2/test1/test-header-1.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=10 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test1/output1_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test2/test-input-2.gra --headerpath=./../testcases/A2/test2/test-header-2.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=8 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test2/output2_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test3/test-input-3.gra --headerpath=./../testcases/A2/test3/test-header-3.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=5 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test3/output3_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test4/test-input-4.gra --headerpath=./../testcases/A2/test4/test-header-4.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=6 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test4/output4_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test5/test-input-5.gra --headerpath=./../testcases/A2/test5/test-header-5.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=8 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test5/output5_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test6/test-input-6.gra --headerpath=./../testcases/A2/test6/test-header-6.dat --outputpath=outputfile.txt --verbose=1 --startk=1 --endk=29 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test6/output6_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test7/test-input-7.gra --headerpath=./../testcases/A2/test7/test-header-7.dat --outputpath=outputfile.txt --verbose=1 --startk=10 --endk=25 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test7/output7_verbose.txt
cd ../COL380-Parallel
time mpirun -np 4 ./a2 --taskid=1 --inputpath=./../testcases/A2/test8/test-input-8.gra --headerpath=./../testcases/A2/test8/test-header-8.dat --outputpath=outputfile.txt --verbose=1 --startk=2 --endk=6 --p=2
cd ../mytests
./a.out ../COL380-Parallel/outputfile.txt ../testcases/A2/test8/output8_verbose.txt
cd ../COL380-Parallel