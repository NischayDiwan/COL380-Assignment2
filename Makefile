all:
	mpic++ main.cpp -std=c++17
clean:
	rm a.out
run:
	mpirun -np 4 ./a.out --taskid=1 --inputpath=./testcases/test3/test-input-3.gra --headerpath=./testcases/test3/test-header-3.dat --outputpath=outputfile.txt --verbose=0 --startk=1 --endk=6 --p=2
