#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);
	int id, sz;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);

	string header = argv[3];
	header = header.substr(13, header.length - 13);	// check 13 ?!
	string filename = argv[2];
	string = string.substr(12, string.length - 12);

	int startk = atoi(argv[6]), endk = atoi(argv[7]);

	// print the above values on sample file	

	MPI_Finalize();
	return 0;
}