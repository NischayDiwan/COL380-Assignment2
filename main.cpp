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
	/*for(auto e: supp){
		pair<int, int> x = e.first;
		cout << x.first << " " << x.second << " " << e.second.size() << endl;
	}*/

	set<pair<int, pair<int, int>>> active;
	map<pair<int, int>, int> hashtable;
	vector<pair<pair<int, int>, int >> T;

	for(auto e: supp){
		pair<int, int> x = e.first;
		active.insert({e.second.size() + 2, x});
		hashtable.insert({x, e.second.size() + 2});
	}

	bool done = 0, action = 0;
	if(active.size() == 0){
		done = 1;
	}
	MPI_Allreduce(&done, &action, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
	action = 1 - action;

	while(action){
		vector<pair<int, int>> cur;
		int min = 100;
		if(active.size()!=0){
			min = (*active.begin()).first;
		}
		int global_min = 100;
		MPI_Allreduce(&min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		if(!done){
			if(active.size() == 0)
				done = 1;
			set<pair<int, pair<int, int>>>::iterator itr1 = active.begin(), itr2 = active.begin();
			for(; itr2 != active.end(); itr2++){
				if((*itr2).first > global_min){
					break;
				}
				else{
					pair<int, int> x = (*itr2).second;
					cur.push_back(x);
					T.push_back({x, (*itr2).first});
					hashtable.erase(x);
				}
			}
			active.erase(itr1, itr2);
			
		}
		vector<pair<pair<int, int>, int>> W;
		for(auto e: cur){
			for(auto w: supp[{e.first, e.second}]){
				W.push_back({e, w});
			}
		}
		int ptr = 0;
		while(true){
			// dst1 for (u, w) edge  and   dst2 for (v, w) edge
			int b[7], u = 0, v = 0, w = 0, dst1 = 0, dst2 = 0;
			if(ptr >= W.size()){
				b[0] = 0;
			}
			else{
				b[0] = 1;
				u = W[ptr].first.first;
				v = W[ptr].first.second;
				w = W[ptr].second;
				if(prio[u] < prio[w]){
					dst1 = u % sz;
					b[1] = u;
					b[2] = w;
				}
				else{
					dst1 = w %sz;
					b[1] = w;
					b[2] = u;
				}
				if(prio[v] < prio[w]){
					dst2 = v % sz;
					b[4] = v;
					b[5] = w;
				}
				else{
					dst2 = w % sz;
					b[4] = w;
					b[5] = v;
				}
				b[3] = dst1; b[6] = dst2;
			}
			MPI_Request msreq;
			MPI_Isend(b, 7, MPI_INT, 0, id, MPI_COMM_WORLD, &msreq);
			int num_packets[sz], payload = 0;
			for(i = 0; i < sz; i++){
				num_packets[i] = 0;
			}
			int mb[7*sz];
			if(id == 0){
				for(i = 0; i < sz; i++){
					MPI_Status status;
					MPI_Recv(mb+7*i, 7, MPI_INT, i, i, MPI_COMM_WORLD, &status);
					if(mb[7*i] != 0){
						num_packets[mb[7*i+3]]++;
						num_packets[mb[7*i+6]]++;
					}
				}
				for(i = 0; i < sz; i++){
					MPI_Request req;
					MPI_Isend(&num_packets[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &req);
				}
			}

			MPI_Status status;
			MPI_Recv(&payload, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &status);

			if(id == 0){
				for(i = 0; i < sz; i++){
					if(mb[7*i] == 1){
						MPI_Request req1, req2;
						MPI_Isend(mb+7*i+1, 2, MPI_INT, mb[7*i+3], mb[7*i+3], MPI_COMM_WORLD, &req1);
						MPI_Isend(mb+7*i+4, 2, MPI_INT, mb[7*i+6], mb[7*i+6], MPI_COMM_WORLD, &req2);
					}
				}
			}

			vector<pair<int, int>> proc;
			int rec[2];
			for(i = 0; i < payload; i++){
				MPI_Status stat;
				MPI_Recv(rec, 2, MPI_INT, 0, id, MPI_COMM_WORLD, &stat);
				proc.push_back({rec[0], rec[1]});
			}
			
			for(auto e: proc){
				if(hashtable.find(e) != hashtable.end()){
					auto erase_itr = active.find({hashtable[e], e});
					active.erase(erase_itr);
					hashtable[e]--;
					active.insert({hashtable[e], e});
				}
			}


			ptr++;
			// break out when all have mb[7k] = 0
			int b_break = 0;
			MPI_Allreduce(&b[0], &b_break, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if(b_break == 0)
				break;
		}

		MPI_Allreduce(&done, &action, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
		action = 1 - action;
	}

	// time measure
	endt = MPI_Wtime();
	if(id == 0){
		cout << "Time taken: " << endt - startt << "\n";
	}
	MPI_Finalize();

	for(auto truss: T){
		pair<int, int> e = truss.first;
		cout << "edge " << e.first << " " << e.second << " has truss number " << e.second;
		cout << endl;
	}
	return 0;
}