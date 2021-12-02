#ifndef __HKPBC_H
#define __HKPBC_H
#include <vector>


const unsigned SIZE = 103968;

class HKPBC
{
public:
	int Lx;
	int Ly;
	int size;
	int* labels;
	int  n_labels = 0; 


public:
	
//	HKPBC(const int a, const int b, std::vector<int> visited(SIZE));
	
	HKPBC(int a, int b, int c);
	
	int uf_find(int x);
	
	void uf_pbc(std::vector<int> &matrix);
	
	int uf_union(int x, int y);

	int uf_make_set(void);

	void uf_initialize(int max_labels); 

	void uf_done(void); 

	std::vector<int> HK(std::vector<int> &matrix);
	
};
#endif
