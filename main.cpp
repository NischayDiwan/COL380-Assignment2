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
	/*int sum = 0;
	for(auto e: supp){
		pair<int, int> x = e.first;
		cout << x.first << " " << x.second << " " << e.second.size() << endl;
		sum += e.second.size();
	}*/

	set<pair<int, pair<int, int>>> active;
	map<pair<int, int>, int> hashtable;
	vector<pair<pair<int, int>, int >> T;
	set<pair<pair<int, int>, int>> Y;

	for(auto e: supp){
		pair<int, int> x = e.first;
		active.insert({e.second.size(), x});
		hashtable.insert({x, e.second.size()});
	}

	int done = 0, action = 0;
	if(active.size() == 0){
		done = 1;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&done, &action, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
	action = 1 - action;
	int loop_cnt = 3;

	while(action!=0){
		vector<pair<int, int>> cur;
		int global_min = loop_cnt - 2;

		/*cout << "NEW ITERTION" << endl;
		cout << "CURRENT GMIN = " << global_min << endl;
		for(auto x: active){
			cout << x.first << " " << x.second.first << " " << x.second.second << endl;
		}
		cout << endl << endl;*/

		if(!done){
			auto itr1 = active.begin(), itr2 = active.begin();
			for(; itr2 != active.end(); itr2++){
				if((*itr2).first >= global_min){
					break;
				}
				else{
					pair<int, int> x = (*itr2).second;
					cur.push_back(x);
					T.push_back({x, loop_cnt - 1});
					hashtable.erase(x);
				}
			}
			active.erase(itr1, itr2);
			if(active.size() == 0)
				done = 1;
		}

		vector<pair<pair<int, int>, int>> W;
		for(auto e: cur){
			for(auto w: supp[{e.first, e.second}]){
				W.push_back({e, w});
			}
		}
		//cout << W.size() << endl;
		/*for(auto x: W){
			cout << "( " << x.first.first << " , " << x.first.second << " , " << x.second << " ) ";
		}*/
		//cout << endl;
		int ptr = 0;
		while(true){
			// dst1 for (u, w) edge  and   dst2 for (v, w) edge
			int b[9], u = 0, v = 0, w = 0, dst1 = 0, dst2 = 0;
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
					b[3] = v;
				}
				else{
					dst1 = w %sz;
					b[1] = w;
					b[2] = u;
					b[3] = v;
				}
				if(prio[v] < prio[w]){
					dst2 = v % sz;
					b[5] = v;
					b[6] = w;
					b[7] = u;
				}
				else{
					dst2 = w % sz;
					b[5] = w;
					b[6] = v;
					b[7] = u;
				}
				b[4] = dst1; b[8] = dst2;
			}
			//parallel
			// MPI_Request msreq;
			// MPI_Isend(b, 9, MPI_INT, 0, id, MPI_COMM_WORLD, &msreq);

			// int num_packets[sz], payload = 0;
			// for(i = 0; i < sz; i++){
			// 	num_packets[i] = 0;
			// }
			// int mb[9*sz];
			// if(id == 0){
			// 	for(i = 0; i < sz; i++){
			// 		MPI_Status status;
			// 		MPI_Recv(mb+9*i, 9, MPI_INT, i, i, MPI_COMM_WORLD, &status);
			// 		if(mb[9*i] != 0){
			// 			num_packets[mb[9*i+4]]++;
			// 			num_packets[mb[9*i+8]]++;
			// 		}
			// 	}
			// 	for(i = 0; i < sz; i++){
			// 		MPI_Request req;
			// 		MPI_Isend(&num_packets[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &req);
			// 	}
			// }

			// MPI_Status status;
			// MPI_Recv(&payload, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &status);

			// gather scatter approach
			int num_packets[sz], payload = 0;
			for(i = 0; i < sz; i++){
				num_packets[i] = 0;
			}
			int mb[9*sz];
			MPI_Gather(b, 9, MPI_INT, mb, 9, MPI_INT, 0, MPI_COMM_WORLD);
			// MPI_Barrier(MPI_COMM_WORLD);
			if(id == 0){
				for(i = 0; i < sz; i++){
					if(mb[9*i] == 1){
						num_packets[mb[9*i+4]]++;
						num_packets[mb[9*i+8]]++;
					}
				}
			}
			MPI_Scatter(num_packets, 1, MPI_INT, &payload, 1, MPI_INT, 0, MPI_COMM_WORLD);

			int tagid[sz];
			for(i = 0; i < sz; i++){
				tagid[i] = 0;
			}
			if(id == 0){
				for(i = 0; i < sz; i++){
					if(mb[9*i] == 1){
						MPI_Request req1, req2;
						MPI_Isend(mb+9*i+1, 3, MPI_INT, mb[9*i+4], mb[9*i+4] + tagid[mb[9*i+4]]++, MPI_COMM_WORLD, &req1);
						MPI_Isend(mb+9*i+5, 3, MPI_INT, mb[9*i+8], mb[9*i+8] + tagid[mb[9*i+8]]++, MPI_COMM_WORLD, &req2);
						// if dst == 0:
						//		do work
					}
				}
			}

			vector<pair<pair<int, int>, int>> proc;
			int rec[3];
			for(i = 0; i < payload; i++){
				MPI_Status stat;
				MPI_Recv(rec, 3, MPI_INT, 0, id + i, MPI_COMM_WORLD, &stat);
				proc.push_back({{rec[0], rec[1]}, rec[2]});
			}

			for(auto x: proc){
				pair<int, int> e = x.first;
				int arr[] = {e.first, e.second, x.second};
				//sort(arr, arr+3);
				if(hashtable.find(e) != hashtable.end()){
					//cout << arr[0] << " " << arr[1] << " " << arr[2] << endl;
					if(Y.find({{arr[0], arr[1]}, arr[2]}) == Y.end()){
						auto erase_itr = active.find({hashtable[e], e});
						active.erase(erase_itr);
						hashtable[e]--;
						active.insert({hashtable[e], e});
						Y.insert({{arr[0], arr[1]}, arr[2]});
						/*if(arr[0] == 1){
							cout << "BOCCHI " << arr[1] << " " << arr[2] << endl;
						}*/
					}
				}
			}

			ptr++;
			// break out when all have mb[7k] = 0
			int b_break = 0;
			//MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&b[0], &b_break, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if(b_break == 0)
				break;
		}

		//MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&done, &action, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
		action = 1 - action;
		// need to map all reduce the below
		int wsize = W.size();
		int decision = 0;
		MPI_Allreduce(&wsize, &decision, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if(decision == 0)
			loop_cnt++;
	}

	// time measure
	endt = MPI_Wtime();
	if(id == 0){
		cout << "Time taken: " << endt - startt << "\n";
	}
	MPI_Finalize();

	// for(auto truss: T){
	// 	pair<int, int> e = truss.first;
	// 	if(truss.second >= 0){
	// 		cout << "edge " << e.first << " " << e.second << " has truss number " << truss.second - 2;
	// 		cout << endl;
	// 	}
	// }
	return 0;
}