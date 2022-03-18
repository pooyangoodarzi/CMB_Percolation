#include"HKPBC_Serial.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <vector>
#include <random>
#include <fstream> 
#include <sstream>
#include <algorithm> 
#include <ctime>  
#include <stdio.h>
#include <stdlib.h>
#include <numeric>

#define maxx(a,b) (a>b?a:b)
#define minn(a,b) (a>b?b:a)

using namespace std;


HKPBC :: HKPBC(int a, int b, int c) 
{
	Lx = a;
	Ly = b;
	size = c;
}


int HKPBC:: uf_find(int x) 
	{
	  int y = x;
	  while (labels[y] != y)
	    y = labels[y];
	
	  while (labels[x] != x) 
	  {
	    int z = labels[x];
	    labels[x] = y;
	    x = z;
	  }
	  return y;
	}
	

int HKPBC:: uf_union(int x, int y) 
	{
	  return labels[uf_find(x)] = uf_find(y);
	}


// int HKPBC:: uf_union_g(int x, int y, std::vector<int> &matrix) 
// 	{
// 		return matrix[uf_find_g(x,matrix)] = uf_find_g(y,matrix);
// 	}


void HKPBC:: union_g(vector<vector<int>> total, int g, int f, int p, int q)
	{
		total[g][p] = uf_union(minn(total[g][p], total[f][q]), maxx(total[g][p], total[f][q]));
		total[f][q] = uf_union(minn(total[g][p], total[f][q]), maxx(total[g][p], total[f][q]));
	}


int HKPBC:: uf_make_set(void) 
	{
	  labels[0] ++;
	  assert(labels[0] < n_labels);
	  labels[labels[0]] = labels[0];
	  return labels[0];
	}
	

void HKPBC:: uf_initialize(int total_size, int max_labels) 
	{
	  n_labels = total_size * max_labels;
	  labels = static_cast<int*>(calloc(sizeof(int), total_size*max_labels));
	  labels[0] = 0;
	}


	
void HKPBC:: uf_done(void) 
	{
	  n_labels = 0;
	  free(labels);
	  labels = 0;
	}
	
vector<int> HKPBC::HK(vector<vector<int>> &matrix, vector<vector<int>> &niegh) //vector<int>
	{
		uf_initialize(matrix.size(), Ly*Lx/2);

		for (int k = 0; k < matrix.size(); k++)
		{
			for (int i=0; i<Lx; i++)
			{
				for (int j=0; j<Ly; j++)
				{	
					if (matrix[k][Ly * i + j]) // if occupied ...
					{                        
						int up = (i==0 ? 0 : matrix[k][Ly * (i - 1) + j]);    //  look up  
						int left = (j==0 ? 0 : matrix[k][Ly * i + (j - 1)]);  //  look left
						
						switch (!!up + !!left)		// direction in which an occupied site exists = 1, else 0
						{
							case 0:
							matrix[k][Ly * i + j] = uf_make_set();      // a new cluster
							break;
							
							case 1:                              // part of an existing cluster
							matrix[k][Ly * i + j] = max(up,left);       // whichever is nonzero is labelled
							break;
							
							case 2:                              // this site binds two clusters
							matrix[k][Ly * i + j] = uf_union(up, left);
							break;
						}
					}
				}
			}
		}
	
	// cout<<"Hk_1:"<<endl;
	// for (int k = 0; k < matrix.size(); k++)
	// {
	// 	cout<<"matrix: "<<k<<endl;
	// 	for(int f=0; f<Lx; f++)
	// 	{
	// 		for(int g=0; g<Ly; g++)
	// 		{
	// 			cout<<matrix[k][f*Ly + g]<<"  ";
	// 		}
	// 		cout<<endl;
	// 	}
	// 	cout<<endl;
	// }
		
	for (int i=0; i<niegh.size(); i++) //niegh.size()
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				if(niegh[i][j] != 100)
				{
					if(niegh[niegh[i][j]][k] == i)
					{
						Hks_pbc(matrix, j, k, i, niegh[i][j]);
					}
				}
			}
		}
	}

	// cout<<"HK_2:"<<endl;
	// for (int k = 0; k < matrix.size(); k++)
	// {
	// 	cout<<"matrix: "<<k<<endl;
	// 	for(int f=0; f<Lx; f++)
	// 	{
	// 		for(int g=0; g<Ly; g++)
	// 		{
	// 			cout<<matrix[k][f*Ly + g]<<"  ";
	// 		}
	// 		cout<<endl;
	// 	}
	// 	cout<<endl;
	// }
	// for(int k = 0; k < matrix.size(); k++)
	// {
	// 	for (int i=0; i<Lx; i++)
	// 	{
	// 		for (int j=0; j<Ly; j++)
	// 		{cout<<labels[matrix[k][Ly * i + j]]<<endl;}
	// 	}
		
	// }

	//   uf_pbc(matrix);   //applying the periodic boundary condition
	
		int *new_labels = static_cast<int*>(calloc(sizeof(int), (matrix.size() * Ly * Lx)/2)); // allocate array, initialized to zero
		new_labels[0] = 0;

		for (int k = 0; k < matrix.size(); k++)
		{
			for (int i=0; i<Lx; i++)
			{
				for (int j=0; j<Ly; j++)
				{ 
					if (matrix[k][Ly * i + j]) 
					{ 
						int x = uf_find(matrix[k][Ly * i + j]);
						if (new_labels[x] == 0) 
						{
							new_labels[0]++;
							new_labels[x] = new_labels[0];
						}
						matrix[k][Ly * i + j] = new_labels[x];
					}
				}
			}
		}
	  


	// cout<<"HK_3:"<<endl;
	// for (int k = 0; k < matrix.size(); k++)
	// {
	// 	cout<<"matrix: "<<k<<endl;
	// 	for(int f=0; f<Lx; f++)
	// 	{
	// 		for(int g=0; g<Ly; g++)
	// 		{
	// 			cout<<matrix[k][f*Ly + g]<<"  ";
	// 		}
	// 		cout<<endl;
	// 	}
	// 	cout<<endl;
	// }
	


	int x = std::accumulate(matrix.begin(), matrix.end(), matrix[0][0],[](int max, const std::vector<int> &v)
								{
									return std::max(max,*std::max_element(v.begin(),v.end()));
								});

	// cout<<"x: "<<x<<endl;
	vector<int> clusters(x + 1, 0);
	vector<int> clus;
	for(int j=0; j<matrix.size(); j++)
	{
		for(int i=0;i<(Lx*Ly);i++)
		{  
			clusters[matrix[j][i]]++; //total[j][i]
		}
	}

	copy(clusters.begin() + 1, clusters.end(), back_inserter(clus));

	// cout<<"the orginal:"<<endl;
	// for(int j=0; j<max_lab + 1; j++)
	// {
	// 	cout<<clusters[j]<<",";
	// }

	// cout<<endl;

	// int nn = sizeof(clus)/sizeof(clus[0]);
	// cout<<"nn: "<<nn<<endl;
	// cout<<"cluster size"<<endl;
	// for (int i = 0; i < nn ; i++)
	// {
	// 	cout<<clus[i]<<",";
	// }
	// cout<<endl;
	
	vector<int> cluster_index(x);
	std::iota(cluster_index.begin(),cluster_index.end(),0); //Initializing
	sort( cluster_index.begin(),cluster_index.end(), [&](int i,int j){return clus[i]<clus[j];} );
	std::for_each(cluster_index.begin(), cluster_index.end(), [](int& d) { d+=1;});

	// cout<<"cluster_index"<<endl;
	// for(int j=0; j<max_lab; j++)
	// {
	// 	cout<<cluster_index[j]<<",";
	// }

	// cout<<endl;

	uf_done();
	return cluster_index;

	}

void HKPBC:: Hks_pbc(vector<vector<int>> &total, int j, int k, int g, int f)
{
	int temp,temp1;
	int n = Lx*Ly;
	if(j==0)
	{
		if(k==0)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<g<<" "<<"B: "<<f<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i] + (g*n)]<<" "<< labels[B[i] + (f*n)]<<endl;
				if(total[g][i] != 0 && total[f][i] != 0 && labels[total[g][i]]!= labels[total[f][i]])
				{
					union_g(total, g, f, i, i);
				}
			}
		}

		if(k==1)
		{
			// cout<<"j: "<<j<<endl; Up
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i] + (g*n)]<<" "<<labels[B[i + (Lx-1)*Ly] + (f*n)]<<endl;
				if(total[g][i] != 0 && total[f][i + (Lx-1)*Ly] != 0 && labels[total[g][i]]!= labels[total[f][i + (Lx-1)*Ly]])
				{
					union_g(total, g, f, i, i + (Lx-1)*Ly);
				}
			}
		}
		if(k==2)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i] + (g*n)]<<" "<<labels[B[(Lx*i)] + (f*n)]<<endl;
				if(total[g][i] != 0 && total[f][(Lx*i)] != 0 && labels[total[g][i]]!= labels[total[f][(Lx*i)]])
				{
					union_g(total, g, f, i, Lx*i);
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i] + (g*n)]<<" "<<labels[B[(Lx*i) + (Ly-1)] + (f*n)]<<endl;
				if(total[g][i] != 0 && total[f][(Lx*i) + (Ly-1)] != 0 && labels[total[g][i]]!= labels[total[f][(Lx*i) + (Ly-1)]])
				{
					union_g(total, g, f, i, (Lx*i) + (Ly-1));
				}
			}
		}
	}

	if(j==1)
	{

		if(k==0)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"i: "<<i<<endl;
				// cout<<"A: "<<g<<" "<<"B: "<<f<<" jk: "<<j<<k<<endl;
				// cout<<A[i + (Lx-1)*Ly]<<" "<<B[i]<<" REAL:"<<A[i + (Lx-1)*Ly] + (g*n)<<" "<<labels[i + (Lx-1)*Ly + (g*n)]<<" REAL: "<<i + (f*n)<<" "<<labels[B[i] + (f*n)]<<endl;
				if(total[g][i + (Lx-1)*Ly] != 0 && total[f][i] != 0 && labels[total[g][i + (Lx-1)*Ly]]!= labels[total[f][i]])
				{
					union_g(total, g, f, i + (Lx-1)*Ly, i);
				}
			}
		}

		if(k==1)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i + (Lx-1)*Ly] + (g*n)]<<" "<<labels[B[i + (Lx-1)*Ly] + (f*n)]<<endl;
				if(total[g][i + (Lx-1)*Ly] != 0 && total[f][i + (Lx-1)*Ly] != 0 && labels[total[g][i + (Lx-1)*Ly]]!= labels[total[f][i + (Lx-1)*Ly]])
				{
					union_g(total, g, f, i + (Lx-1)*Ly, i + (Lx-1)*Ly);
				}
			}
		}
		if(k==2)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i + (Lx-1)*Ly] + (g*n)]<<" "<<labels[B[(Lx*i)] + (f*n)]<<endl;
				if(total[g][i + (Lx-1)*Ly] != 0 && total[f][(Lx*i)] != 0 && labels[total[g][i + (Lx-1)*Ly]]!= labels[total[f][(Lx*i)]])
				{
					union_g(total, g, f, i + (Lx-1)*Ly, (Lx*i));
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[i + (Lx-1)*Ly] + (g*n)]<<" "<<labels[B[(Lx*i) + (Ly-1)] + (f*n)]<<endl;
				if(total[g][i + (Lx-1)*Ly] != 0 && total[f][(Lx*i) + (Ly-1)] != 0 && labels[total[g][i + (Lx-1)*Ly]]!= labels[total[f][(Lx*i) + (Ly-1)]])
				{
					union_g(total, g, f, i + (Lx-1)*Ly, (Lx*i) + (Ly-1));
				}
			}
		}
	}

	if(j==2)
	{
		if(k==0)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[Lx*i] + (g*n)]<<" "<<labels[B[i] + (f*n)]<<endl;
				if(total[g][Lx*i] != 0 && total[f][i] != 0 && labels[total[g][Lx*i]]!= labels[total[f][i]])
				{
					union_g(total, g, f, Lx*i, i);
				}
			}
		}

		if(k==1)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[Lx*i] + (g*n)]<<" "<<labels[B[i + (Lx-1)*Ly] + (f*n)]<<endl;
				if(total[g][Lx*i] != 0 && total[f][i + (Lx-1)*Ly] != 0 && labels[total[g][Lx*i]]!= labels[total[f][i + (Lx-1)*Ly]])
				{
					union_g(total, g, f, Lx*i, i + (Lx-1)*Ly);
				}
			}
		}
		if(k==2)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[Lx*i] + (g*n)]<<" "<<labels[B[(Lx*i)] + (f*n)]<<endl;
				if(total[g][Lx*i] != 0 && total[f][Lx*i] != 0 && labels[total[g][Lx*i]]!= labels[total[f][Lx*i]])
				{
					union_g(total, g, f, Lx*i, Lx*i);
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<g<<" "<<"B: "<<f<<" jk: "<<j<<k<<endl;
				// cout<<labels[Lx*i + (g*n)]<<" "<<labels[(Lx*i) + (Ly-1) + (f*n)]<<endl;	
				if(total[g][Lx*i] != 0 && total[f][(Lx*i) + (Ly-1)] != 0 && labels[total[g][Lx*i]] != labels[total[f][(Lx*i) + (Ly-1)]])
				{
					union_g(total, g, f, Lx*i, (Lx*i) + (Ly-1));
				}
			}
		}
	}

	if(j==3)
	{
		if(k==0)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[(Lx*i) + (Ly-1)] + (g*n)]<<" "<<labels[B[i] + (f*n)]<<endl;
				if(total[g][(Lx*i) + (Ly-1)] != 0 && total[f][i] != 0 && labels[total[g][(Lx*i) + (Ly-1)]]!= labels[total[f][i]])
				{
					union_g(total, g, f, (Lx*i) + (Ly-1), i);
				}
			}
		}

		if(k==1)
		{
			for(int i=0; i<Ly; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[(Lx*i) + (Ly-1)] + (g*n)]<<" "<<labels[B[i + (Lx-1)*Ly] + (f*n)]<<endl;
				if(total[g][(Lx*i) + (Ly-1)] != 0 && total[f][i + (Lx-1)*Ly] != 0 && labels[total[g][(Lx*i) + (Ly-1)]]!= labels[total[f][i + (Lx-1)*Ly]])
				{
					union_g(total, g, f, (Lx*i) + (Ly-1), i + (Lx-1)*Ly);
				}
			}
		}
		if(k==2)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[(Lx*i) + (Ly-1)] + (g*n)]<<" "<<labels[B[(Lx*i)] + (f*n)]<<endl;
				if(total[g][(Lx*i) + (Ly-1)] != 0 && total[f][(Lx*i)] != 0 && labels[total[g][(Lx*i) + (Ly-1)]]!= labels[total[f][(Lx*i)]])
				{
					union_g(total, g, f, (Lx*i) + (Ly-1), Lx*i);
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				// cout<<"A: "<<f<<" "<<"B: "<<g<<" jk: "<<j<<k<<endl;
				// cout<<labels[A[(Lx*i) + (Ly-1)] + (g*n)]<<" "<<labels[B[(Lx*i) + (Ly-1)] + (f*n)]<<endl;
				if(total[g][(Lx*i) + (Ly-1)] != 0 && total[f][(Lx*i) + (Ly-1)] != 0 && labels[total[g][(Lx*i) + (Ly-1)]]!= labels[total[f][(Lx*i) + (Ly-1)]])
				{
					union_g(total, g, f, (Lx*i) + (Ly-1), (Lx*i) + (Ly-1));
				}
			}
		}
	}

} 

