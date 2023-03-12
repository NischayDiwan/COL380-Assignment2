#include<bits/stdc++.h>
//#include<mpi.h>

using namespace std;

// reads value from file into parameters
void read(string filename, int &n, int &m, vector<int> &val, vector<int> &A){
	ifstream inFile;
	inFile.open(filename, ios::binary | ios::in);

	inFile.read(reinterpret_cast<char *>(&n), 4);
	inFile.read(reinterpret_cast<char *>(&m), 4);

	for(int i = 0; i < n; i++){
		int deg = 0, tmp;
		inFile.read(reinterpret_cast<char *>(&tmp), 4);
		inFile.read(reinterpret_cast<char *>(&deg), 4);
		val.push_back(tmp);
		for(int j = 0; j < deg; j++){
			inFile.read(reinterpret_cast<char *>(&tmp), 4);
			A.push_back(tmp);
		}
	}
	inFile.close();
}

int main(int argc, char* argv[]){
	int n, m, i, j;
	vector<int> val, A;

	read("", n, m, val, A);

	vector<int> adj[n];

	return 0;
}