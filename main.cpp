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
	int startk = stoi(stk), endk = stoi(enk);

	// store number of nodes and edges
	ifstream infile(filename, ios::in | ios::binary);
	infile.read(reinterpret_cast<char *>(&n), 4);
	infile.read(reinterpret_cast<char *>(&m), 4);
	infile.close();

	pair<int, int> deg[n]; // store ( deg[v], v )
	int offset[n], prio[n];	// stores the offset of i'th node

	// vertex i starts at 4*i bytes in header
	ifstream hfile(header, ios::in | ios::binary);
	for(i = 0; i < n; i++){
		hfile.read(reinterpret_cast<char *>(offset + i), 4);	// stores offset
	}
	hfile.close();

	infile.open(filename, ios::in | ios::binary);
	for(i = 0; i < n; i++){
		infile.seekg(offset[i] + 4, ios::beg);
		infile.read(reinterpret_cast<char *>(&tmp), 4);
		deg[i] = {tmp, i};
	}
	infile.close();

	// establish order on vertices based on degree
	sort(deg, deg + n);
	for(i = 0; i < n; i++){
		prio[deg[i].second] = n - i;
		//cout << deg[i].first << " " << deg[i].second << "\n";
	}

	// assign vertices to the nodes [ based modulo size ]
	// node i gets vertices [i, sz + i, 2*sz + i, ...]  ----> check stupid cases ( sz = 1 , seems ok)
	// edge e = (u, v) then e belongs to processor having lower priority vertex
	vector<int> V;
	for(i = id; i < n; i += sz){
		V.push_back(i);
	}

	int num_nodes = V.size();
	vector<int> E[num_nodes];	// E[u] will stores edges (u, v) with u < v

	infile.open(filename, ios:: in | ios::binary);
	for(auto v: V){
		infile.seekg(offset[v] + 4, ios::beg);
		infile.read(reinterpret_cast<char *>(&tmp), 4);
		for(i = 0; i < tmp; i++){
			infile.read(reinterpret_cast<char *>(&j), 4);
			if(prio[j] > prio[v]){
				// cout << v << " " << j << endl;
				E[(v - id)/sz].push_back(j);
			}
		}
	}
	infile.close();

	// cout << id << " " << sz << endl;
	for(i = 0; i < num_nodes; i++){
		string output = "";
		for(auto u: E[i]){
			cout << u << " ";
			output += u;
		}
		cout << endl;
	}

	// triangle enumeration



	MPI_Finalize();
	return 0;
}