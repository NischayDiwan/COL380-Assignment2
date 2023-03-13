all: main.cpp
	mpic++ main.cpp -std=c++17 -o a2
run: all
	mpirun -np 2 ./a2 --taskid=1 --inputpath=./testcases/test4/test-input-4.gra --headerpath=./testcases/test4/test-header-4.dat --outputpath=outputfile.txt --verbose=0 --startk=1 --endk=6 --p=2
clean:
	rm a2