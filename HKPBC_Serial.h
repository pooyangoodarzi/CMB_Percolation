#ifndef __HKPBC_H
#define __HKPBC_H
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
	
	void uf_pbc(vector<int> &matrix);
	
	int uf_union(int x, int y);

	// int uf_union_g(int x, int y, vector<int> &matrix); 

	int uf_make_set(void);

	void uf_initialize(int max_labels, int last_label, int rank);

	void uf_done(void); 

	int HK(vector<int> &matrix, int rank, int last_label);

	vector<int> serial(vector<vector<int>> &niegh, vector<vector<int>> &total, int max_lab);

	void uf_Ginitialize(int max_labels_g);

	void Hks_pbc(vector<int> &A, vector<int> &B, int j);

	// int uf_find_g(int x, std::vector<int> &matrix);

	
};
#endif
