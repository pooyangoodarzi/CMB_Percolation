#ifndef __HKPBC_Serial_H
#define __HKPBC_Serial_H
#include <vector>

using namespace std;

const unsigned SIZE = 103968;

class HKPBC
{
public:
	int Lx;
	int Ly;
	int num;
	int* labels;
	int* labels_g;
	int  n_labels = 0; 
	int  n_labels_g = 0; 
	int size, rank, last_label;


public:
	
//	HKPBC(const int a, const int b, std::vector<int> visited(SIZE));
	
	HKPBC(int a, int b, int c);
	
	int uf_find(int x);
	
	// void uf_pbc(vector<int> &matrix);
	
	int uf_union(int x, int y);

	int uf_make_set(void); //void

	void uf_initialize(int total_size, int max_labels);

	void uf_done(void); 

	vector<int> HK(vector<vector<int>> &matrix, vector<vector<int>> &niegh);

	void Hks_pbc(vector<vector<int>> &total, int j, int k, int g, int f);

	void union_g(vector<vector<int>> total, int g, int f, int p, int q);
	
};
#endif
