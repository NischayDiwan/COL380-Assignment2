#include<bits/stdc++.h>
#include<mpi.h>

using namespace std;

int main(int argc, char* argv[]){

	int i, j, k, n, m, tmp;
	// time measure variables
	double startt, endt;

	MPI_Init(&argc, &argv);
	int id, sz;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);
	startt = MPI_Wtime();

	// read input arguments
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
		prio[deg[i].second] = i;
		//cout << deg[i].first << " " << deg[i].second << "\n";
		//cout << deg[i].second << " ";
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
	vector<int> E[num_nodes];	// E[u] will stores edges (u, v) 
	unordered_set<int> target;
	map<int, vector<int> > par;
	infile.open(filename, ios:: in | ios::binary);
	int count = 0;
	for(auto v: V){
		infile.seekg(offset[v] + 4, ios::beg);
		infile.read(reinterpret_cast<char *>(&tmp), 4);
		for(i = 0; i < tmp; i++){
			infile.read(reinterpret_cast<char *>(&j), 4);
			E[count].push_back(j);
			if(prio[j] > prio[v]){
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
	map<pair<int, int>, vector<int>> supp;

	// triangle enumeration
	infile.open(filename, ios::in | ios::binary);
	// iterate over (u, v, w) with u < v and u < w. Check if vw exists
	// first loop over all v? so that we can read into vector --> then loop over vertices of current processor? 
	map<pair<pair<int, int>, int>, bool> processed;
	for(auto v: target){
		vector<int> adj;
		infile.seekg(offset[v] + 4, ios::beg);
		infile.read(reinterpret_cast<char *>(&count), 4);
		for(i = 0; i < count; i++){
			infile.read(reinterpret_cast<char *>(&tmp), 4);
			adj.push_back(tmp);
		}
		for(auto u: par[v]){
			/*if(v <= u)
				continue; */
			for(auto w: E[(u - id)/sz]){
				if(v == w)
					continue;
				else if(binary_search(adj.begin(), adj.end(), w)){
					// increment support
					if(supp.find({u, v})!=supp.end()){
						supp[{u, v}].push_back(w);
					}
					else{
						vector<int> throwaway;
						throwaway.push_back(w);
						supp.insert({{u, v}, throwaway});
					}
				}
			}
		}
	}
	infile.close();
	// supp[u, v] stores support of edge (u, v) in current processor
	//cout << endl << "done counting " << endl << endl;
	for(auto e: supp){
		pair<int, int> x = e.first;
		cout << x.first << " " << x.second << " " << e.second.size() << endl;
	}

	set<pair<int, pair<int, int>>> active;
	vector<pair<pair<int, int>, int >> T;

	for(auto e: supp){
		pair<int, int> x = e.first;
		active.insert({e.second.size() + 2, {x.first, x.second}});
	}

	bool done = false, action = false;
	if(active.size() == 0){
		done = true;
	}

	MPI_Allreduce(&done, &action, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD);
	action = not action;
	
	while(action){
		if(!done){
			int min = (*active.begin()).first;
			set<pair<int, pair<int, int>>>::iterator itrs1 = active.begin(), itr2 = active.begin();
			vector<pair<int, int>> cur;
			for(; itr2 != active.end(); itr2++){
				if((*itr2).first > min)
					break;
				else{
					pair<int, int> x = (*itr2).second;
					cur.push_back(x);
					T.push_back({x, (*itr2).first});
				}
			}
			active.erase(itr1, itr2);

			for(auto e: cur){
				for(auto w: supp[{e.first, e.second}]){
					// send msg w to rank : (w % sz)
					/*int hj;
					MPI_Request req;
					MPI_Isend(&w, 1, MPI_INT, w%sz, id, MPI_COMM_WORLD, &req);
					MPI_Irecv(&hj, 1, MPI_INT, )*/
				}
			}



			if(active.size() == 0)
				done = true;
		}

		MPI_Allreduce(&done, &action, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD);
		action = not action;
	}

	// time measure
	endt = MPI_Wtime();
	if(id == 0){
		cout << "Time taken: " << endt - startt << "\n";
	}
	MPI_Finalize();
	return 0;
}