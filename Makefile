all:
	mpic++ main.cpp
clean:
	rm a.out
run:
	mpirun -np 1 ./a.out --taskid=1 --inputpath=./testcases/test4/test-input-4.gra --headerpath=./testcases/test4/test-header-4.dat --outputpath=outputfile.txt --verbose=0 --startk=1 --endk=6 --p=2