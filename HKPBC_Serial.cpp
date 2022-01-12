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

void HKPBC:: pbc_union(vector<vector<int>> &total, int oldlabel, int newlabel)
{
	for (int i = 0; i < total.size(); i++)
	{
		replace(total[i].begin(), total[i].end(), oldlabel, newlabel);
	}
	
}

	// int HKPBC:: uf_find_g(int x, std::vector<int> &matrix) 
	// {
	// 	cout<<x<<endl;
	// 	int y = x;
	// 	cout<<matrix[y]<<endl;
	// 	while (matrix[y] != y)
	// 	y = matrix[y];
	// 	cout<<"hi1"<<endl;
	// 	while (matrix[x] != x) 
	// 	{
	// 	int z = matrix[x];
	// 	matrix[x] = y;
	// 	x = z;
	// 	}
	// 	cout<<"hi2"<<endl;
	// 	return y;
	// }
	
	
// void HKPBC:: uf_pbc(vector<int> &matrix)
// 	{
// 		for(int i=0; i<Lx; i++)
// 	  {
// 		if(matrix[Ly * i] && matrix[Ly * i + (Ly-1)] && labels[matrix[Ly * i]]!= labels[matrix[Ly * i + (Ly-1)]])
// 		{   
// 			matrix[Ly * i] = uf_union(matrix[Ly * i], matrix[Ly * i + (Ly-1)]);
// 		}
// 	  } 	
// 	}
	

int HKPBC:: uf_union(int x, int y) 
	{
	  return labels[uf_find(x)] = uf_find(y);
	}


// int HKPBC:: uf_union_g(int x, int y, std::vector<int> &matrix) 
// 	{
// 		return matrix[uf_find_g(x,matrix)] = uf_find_g(y,matrix);
// 	}


int HKPBC:: uf_make_set(void) 
	{
	  labels[0] ++;
	  assert(labels[0] < n_labels);
	  labels[labels[0]] = labels[0];
	  return labels[0];
	}
	

void HKPBC:: uf_initialize(int max_labels) 
	{	
	  n_labels = max_labels;
	  labels = static_cast<int*>(calloc(sizeof(int), max_labels));
	  labels[0] = 0;
	}


	
void HKPBC:: uf_done(void) 
	{
	  n_labels = 0;
	  free(labels);
	  labels = 0;
	}
	
int HKPBC::HK(vector<int> &matrix, int last_label) //vector<int>
	{
		uf_initialize((Ly * Lx)/2);
		for (int i=0; i<Lx; i++)
		{
			for (int j=0; j<Ly; j++)
			{	
				if (matrix[Ly * i + j]) // if occupied ...
				{                        
					int up = (i==0 ? 0 : matrix[Ly * (i - 1) + j]);    //  look up  
					int left = (j==0 ? 0 : matrix[Ly * i + (j - 1)]);  //  look left
					
					switch (!!up + !!left)		// direction in which an occupied site exists = 1, else 0
					{
						case 0:
						matrix[Ly * i + j] = uf_make_set();      // a new cluster
						break;
						
						case 1:                              // part of an existing cluster
						matrix[Ly * i + j] = max(up,left);       // whichever is nonzero is labelled
						break;
						
						case 2:                              // this site binds two clusters
						matrix[Ly * i + j] = uf_union(up, left);
						break;
					}
				}
			}
		}
	  

	//   uf_pbc(matrix);   //applying the periodic boundary condition
	
	  int *new_labels = static_cast<int*>(calloc(sizeof(int), (Ly * Lx)/2)); // allocate array, initialized to zero
	  new_labels[0] = last_label;
	  for (int i=0; i<Lx; i++)
	  {
	    for (int j=0; j<Ly; j++)
		{ 
	      if (matrix[Ly * i + j]) 
		  { 
			int x = uf_find(matrix[Ly * i + j]);
			if (new_labels[x] == 0) 
			{
				new_labels[0]++;
				new_labels[x] = new_labels[0];
			}
			matrix[Ly * i + j] = new_labels[x];
	      }
		}
	  }

	  int total_clusters = new_labels[0];	
	//   cout<<"total: "<<total_clusters<<endl;
	  free(new_labels);
	  uf_done();
	  
	//    cout<<"HK:"<<endl;
	//    for(int f=0; f<Lx; f++)
	//  	{
	//  		for(int g=0; g<Ly; g++)
	//  		{
	//  			cout<<matrix[f*Ly + g]<<"  ";
	//  		}
	//  		cout<<endl;
	//  	}
	//    cout<<endl;
	  
	  vector<int> clusters(total_clusters + 1, 0);
	  vector<int> clus;
	  for(int i=0;i<(Lx*Ly);i++)
		{  
			clusters[matrix[i]]++;
		}
	  copy(clusters.begin() + 1, clusters.end(), back_inserter(clus));
	//   cout<<"cluster size"<<endl;
	//   for (int i = 0; i < total_clusters ; i++)
	//   {
	// 	cout<<clus[i]<<",";
	//   }
	//   cout<<endl;
	
		vector<int> cluster_index(total_clusters);
		std::iota(cluster_index.begin(),cluster_index.end(),0); //Initializing
		sort( cluster_index.begin(),cluster_index.end(), [&](int i,int j){return clus[i]<clus[j];} );
		std::for_each(cluster_index.begin(), cluster_index.end(), [](int& d) { d+=1;});
	//   int bigclus = *max_element(clusters.begin(), clusters.end());
	//  cout<<"big cluster: "<<bigclus<<endl;
	//   int data = distance(clusters.begin(), find(clusters.begin(), clusters.end(), bigclus))+1;
	//  cout<<"index:"<<data<<endl;
	//   vector <int> index;
	  
		// cout<<"indexes"<<endl;
		// for(int i=0; i<total_clusters; i++)
		// {
		// 	cout<<cluster_index[i]<<",";
		// }
		// cout<<endl;

		// cout<<"clusters: "<<endl;
		// for(int i=0; i<total_clusters; i++)
		// {
		// 	cout<<clusters[i]<<",";
		// }
		// cout<<endl;

	//   return cluster_index;
	  return total_clusters;
	}

void HKPBC:: Hks_pbc(vector<int> &A, vector<int> &B, vector<vector<int>> &total, int j, int k)
{
	int temp,temp1;
	if(j==0)
	{
		if(k==0)
		{
			for(int i=0; i<Ly; i++)
			{
				if(A[i] != 0 && B[i] != 0 && A[i]!= B[i])
				{
					temp = B[i];
					temp1 = A[i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}

		if(k==1)
		{
			// cout<<"j: "<<j<<endl; Up
			for(int i=0; i<Ly; i++)
			{
				// cout<<i<<","<<A[i]<<" "<<i + (Lx-1)*Ly<<","<<B[i + (Lx-1)*Ly]<<endl;
				if(A[i] != 0 && B[i + (Lx-1)*Ly] != 0 && A[i]!= B[i + (Lx-1)*Ly])
				{
					temp = B[i + (Lx-1)*Ly];
					temp1 = A[i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);

				}
			}
		}
		if(k==2)
		{
			// cout<<"j: "<<j<<endl; Up
			for(int i=0; i<Lx; i++)
			{
				// cout<<i<<","<<A[i]<<" "<<i + (Lx-1)*Ly<<","<<B[i + (Lx-1)*Ly]<<endl;
				if(A[i] != 0 && B[(Lx*i)] != 0 && A[i]!= B[(Lx*i)])
				{
					temp = B[(Lx*i)];
					temp1 = A[i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
				
		if(k==3)
		{
			// cout<<"j: "<<j<<endl; Up
			for(int i=0; i<Lx; i++)
			{
				// cout<<i<<","<<A[i]<<" "<<i + (Lx-1)*Ly<<","<<B[i + (Lx-1)*Ly]<<endl;
				if(A[i] != 0 && B[(Lx*i) + (Ly-1)] != 0 && A[i]!= B[(Lx*i) + (Ly-1)])
				{
					temp = B[(Lx*i) + (Ly-1)];
					temp1 = A[i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
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
				if(A[i + (Lx-1)*Ly] != 0 && B[i] != 0 && A[i + (Lx-1)*Ly]!= B[i])
				{
					temp = B[i];
					temp1 = A[i + (Lx-1)*Ly];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}

		if(k==1)
		{
			for(int i=0; i<Ly; i++)
			{
				if(A[i + (Lx-1)*Ly] != 0 && B[i + (Lx-1)*Ly] != 0 && A[i + (Lx-1)*Ly]!= B[i + (Lx-1)*Ly])
				{
					temp = B[i + (Lx-1)*Ly];
					temp1 = A[i + (Lx-1)*Ly];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
		if(k==2)
		{
			// cout<<"j: "<<j<<endl; Up
			for(int i=0; i<Lx; i++)
			{
				// cout<<i<<","<<A[i]<<" "<<i + (Lx-1)*Ly<<","<<B[i + (Lx-1)*Ly]<<endl;
				if(A[i + (Lx-1)*Ly] != 0 && B[(Lx*i)] != 0 && A[i + (Lx-1)*Ly]!= B[(Lx*i)])
				{
					temp = B[(Lx*i)];
					temp1 = A[i + (Lx-1)*Ly];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				if(A[i + (Lx-1)*Ly] != 0 && B[(Lx*i) + (Ly-1)] != 0 && A[i + (Lx-1)*Ly]!= B[(Lx*i) + (Ly-1)])
				{
					temp = B[(Lx*i) + (Ly-1)];
					temp1 = A[i + (Lx-1)*Ly];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
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
				if(A[Lx*i] != 0 && B[i] != 0 && A[Lx*i]!= B[i])
				{
					temp = B[i];
					temp1 = A[Lx*i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}

		if(k==1)
		{
			for(int i=0; i<Ly; i++)
			{
				if(A[Lx*i] != 0 && B[i + (Lx-1)*Ly] != 0 && A[Lx*i]!= B[i + (Lx-1)*Ly])
				{
					temp = B[i + (Lx-1)*Ly];
					temp1 = A[Lx*i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
		if(k==2)
		{
			for(int i=0; i<Lx; i++)
			{
				if(A[Lx*i] != 0 && B[(Lx*i)] != 0 && A[Lx*i]!= B[(Lx*i)])
				{
					temp = B[(Lx*i)];
					temp1 = A[Lx*i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				if(A[Lx*i] != 0 && B[(Lx*i) + (Ly-1)] != 0 && A[Lx*i]!= B[(Lx*i) + (Ly-1)])
				{
					temp = B[(Lx*i) + (Ly-1)];
					temp1 = A[Lx*i];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
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
				if(A[(Lx*i) + (Ly-1)] != 0 && B[i] != 0 && A[(Lx*i) + (Ly-1)]!= B[i])
				{
					temp = B[i];
					temp1 = A[(Lx*i) + (Ly-1)];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}

		if(k==1)
		{
			for(int i=0; i<Ly; i++)
			{
				if(A[(Lx*i) + (Ly-1)] != 0 && B[i + (Lx-1)*Ly] != 0 && A[(Lx*i) + (Ly-1)]!= B[i + (Lx-1)*Ly])
				{
					temp = B[i + (Lx-1)*Ly];
					temp1 = A[(Lx*i) + (Ly-1)];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
		if(k==2)
		{
			for(int i=0; i<Lx; i++)
			{
				if(A[(Lx*i) + (Ly-1)] != 0 && B[(Lx*i)] != 0 && A[(Lx*i) + (Ly-1)]!= B[(Lx*i)])
				{
					temp = B[(Lx*i)];
					temp1 = A[(Lx*i) + (Ly-1)];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
				
		if(k==3)
		{
			for(int i=0; i<Lx; i++)
			{
				if(A[(Lx*i) + (Ly-1)] != 0 && B[(Lx*i) + (Ly-1)] != 0 && A[(Lx*i) + (Ly-1)]!= B[(Lx*i) + (Ly-1)])
				{
					temp = B[(Lx*i) + (Ly-1)];
					temp1 = A[(Lx*i) + (Ly-1)];
					vector<int> tempo{temp1, temp};
					sort(tempo.begin(), tempo.end(), [](int a, int b) {return a < b;});
					pbc_union(total, tempo[1], tempo[0]);
				}
			}
		}
	}

} 

vector<int> HKPBC:: serial(vector<vector<int>> &niegh, vector<vector<int>> &total, int max_lab)
{
	for (int i=0; i<niegh.size(); i++) //niegh.size()
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				if(niegh[niegh[i][j]][k] == i)
				{
					Hks_pbc(total[i], total[niegh[i][j]], total, j, k);
				}
			}
		}
	}

	// for (int i=0; i<niegh.size(); i++) //niegh.size()
	// {
	// 	for (int j = 0; j < 4; j++)
	// 	{
	// 		// cout<<i<<","<<niegh[i][j]<<endl;
	// 		for (int k = 0; k < 4; k++)
	// 		{
	// 			// cout<<i<<","<<niegh[niegh[i][j]][k]<<","<<k<<","<<endl;
	// 			if(niegh[niegh[i][j]][k] == i)
	// 			{
	// 				Hks_pbc(total[i], total[niegh[i][j]], j, k);
	// 			}
	// 		}
	// 	}
	// }


	int x = std::accumulate(total.begin(), total.end(), total[0][0],[](int max, const std::vector<int> &v)
								{
									return std::max(max,*std::max_element(v.begin(),v.end()));
								});

	// cout<<"x: "<<x<<endl;
	int clusters[x+1] = {0};
	int clus[x];
	for(int j=0; j<total.size(); j++)
	{
		for(int i=0;i<(Lx*Ly);i++)
		{  
			clusters[total[j][i]]++; //total[j][i]
		}
	}

	// cout<<"the orginal:"<<endl;
	// for(int j=0; j<max_lab + 1; j++)
	// {
	// 	cout<<clusters[j]<<",";
	// }

	// cout<<endl;

	copy(clusters + 1, clusters + x +1, clus);

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


