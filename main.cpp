#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[]){

	int i, j, k, n, m, tmp;

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
	// BUT THIS HAS TO BE MADE ACCORDING TO ASCENDING ORDER OF DEGREE
	// edge e = (u, v) then e belongs to processor having lower priority vertex
	vector<int> V;
	for(i = id; i < n; i += sz){
		V.push_back(i);
	}

	int num_nodes = V.size();
	vector<int> E[num_nodes];	// E[u] will stores edges (u, v) with u < v 
	unordered_set<int> target;
	map<int, vector<int> > par;
	infile.open(filename, ios:: in | ios::binary);
	int count = 0;
	for(auto v: V){
		infile.seekg(offset[v] + 4, ios::beg);
		infile.read(reinterpret_cast<char *>(&tmp), 4);
		for(i = 0; i < tmp; i++){
			infile.read(reinterpret_cast<char *>(&j), 4);
			if(prio[j] > prio[v]){
				E[count].push_back(j);
				target.insert(j);
				if(par.find(j)!=par.end()){
					par[j].push_back(v);
				}
				else{
					vector<int> throwaway;
					throwaway.push_back(v);
					par.insert({j, throwaway});
				}
			}
		}
		count += 1;
	}
	infile.close();

	// count of (x, y) in a triangle
	map<pair<int, int>, int> supp;

	// triangle enumeration
	infile.open(filename, ios::in | ios::binary);
	// iterate over (u, v, w) with u < v and u < w. Check if vw exists
	// first loop over all v? so that we can read into vector --> then loop over vertices of current processor? 
	for(auto& v: target){
		vector<int> adj;
		infile.seekg(offset[v] + 4, ios::beg);
		infile.read(reinterpret_cast<char *>(&count), 4);
		for(i = 0; i < tmp; i++){
			infile.read(reinterpret_cast<char *>(&tmp), 4);
			adj.push_back(tmp);
		}
		for(auto u: par[v]){
			for(auto w: E[(u - id)/sz]){
				if(prio[w] < prio[v])
					continue;
				// check if vw is an edge in G
				if(binary_search(adj.begin(), adj.end(), w)){
					// increment support
					if(supp.find({u, v})!=supp.end()){
						supp[{u, v}]++;
					}
					else{
						supp.insert({{u, v}, 1});
					}
					if(supp.find({u, w})!=supp.end()){
						supp[{u, w}]++;
					}
					else{
						supp.insert({{u, w}, 1});
					}
					if(supp.find({v, w})!=supp.end()){
						supp[{v, w}]++;
					}
					else{
						supp.insert({{v, w}, 1});
					}
				}
			}
		}
	}
	infile.close();

	// have to add up the supports present in some other processor


	MPI_Finalize();
	return 0;
}