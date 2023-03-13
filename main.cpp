#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[]){

	int i, j, n, m, tmp;

	MPI_Init(&argc, &argv);
	int id, sz;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);

	string header = argv[3];
	header = header.substr(13, header.length() - 13);	// check 13 ?!
	string filename = argv[2];
	filename = filename.substr(12, filename.length() - 12);

	string stk = argv[6], enk = argv[7];
	stk = stk.substr(9, stk.length() - 9);
	enk = enk.substr(7, enk.length() - 7);
	int startk = atoi(stk), endk = atoi(enk);

	// store number of nodes and edges
	ifstream infile(filename, ios::in | ios::binary);
	infile.read(reinterpret_cast<char *>(&n), 4);
	infile.read(reinterpret_cast<char *>(&m), 4);
	infile.close();

	pair<int, int> deg[n]; // store ( deg[v], v )
	int offset[n];	// stores the offset of i'th node

	// vertex i starts at 4*i bytes in header
	ifstream hfile(header, ios::in | ios::binary);
	for(i = 0; i < n; i++){
		hfile.read(reinterpret_cast<char *>(offset + i), 4);	// stores offset
	}
	hfile.close();

	ifstream infile(filename, ios::in | ios::binary);
	for(i = 0; i < n; i++){
		infile.seekp(offset[i], ios::beg + 4);
		//infile.read(reinterpret_cast<char *>(&tmp), 4);
		infile.read(reinterpret_cast<char *>(&tmp), 4);
		deg[i] = {tmp, i};
	}
	infile.close();

	sort(deg, deg + n);
	for(i = 0; i < n; i++){
		cout << deg[i].first << " " << deg[i].second << "\n";
	}

	MPI_Finalize();
	return 0;
}