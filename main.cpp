#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[]){

	MPI_Init(&argc, &argv);
	int id, sz;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);

	string header = argv[3].substr(13, argv[3].length - 13);	// check 13 ?!
	string filename = argv[2].substr(12, argv[2].length - 12);

	startk = atoi(argv[6]);
	endk = atoi(argv[7]);

	// print the above values on sample file	

	MPI_Finalize();
	return 0;
}