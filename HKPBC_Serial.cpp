#include"HKPBC.h"
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

using namespace std;


HKPBC :: HKPBC(int a, int b, int c) 
{
	Lx = a;
	Ly = b;
	size = c;
	uf_Ginitialize(size / 2);

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
	
	
void HKPBC:: uf_pbc(vector<int> &matrix)
	{
		for(int i=0; i<Lx; i++)
	  {
		if(matrix[Ly * i] && matrix[Ly * i + (Ly-1)] && labels[matrix[Ly * i]]!= labels[matrix[Ly * i + (Ly-1)]])
		{   
			matrix[Ly * i] = uf_union(matrix[Ly * i], matrix[Ly * i + (Ly-1)]);
		}
	  } 	
	}
	

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
	

void HKPBC:: uf_initialize(int max_labels, int last_label, int rank) 
	{	
	  n_labels = (rank * max_labels);
	  labels = static_cast<int*>(calloc(sizeof(int), n_labels));
	  labels[0] = last_label;
	//   cout<<labels[0]<<endl;
	}

void HKPBC:: uf_Ginitialize(int max_labels_g) 
	{	
	  n_labels_g = max_labels_g;
	  labels_g = static_cast<int*>(calloc(sizeof(int), n_labels_g));
	  labels_g[0] = 0;
	}
	
void HKPBC:: uf_done(void) 
	{
	  n_labels = 0;
	  free(labels);
	  labels = 0;
	}
	
int HKPBC::HK(vector<int> &matrix, int rank, int last_label) //vector<int>
	{
	
		uf_initialize((Ly * Lx / 2), last_label, rank);
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
	
	  int *new_labels = static_cast<int*>(calloc(sizeof(int), n_labels)); // allocate array, initialized to zero
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
	  
	  int clusters[total_clusters + 1] = {0};
	  int clus[total_clusters];
	  for(int i=0;i<(Lx*Ly);i++)
		{  
			clusters[matrix[i]]++;
		}
	  copy(clusters + 1, clusters + total_clusters +1, clus);
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

void HKPBC:: Hks_pbc(vector<int> &A, vector<int> &B, int j)
{
	int temp,temp1;
	if(j==0)
	{
		// cout<<"j: "<<j<<endl;
		for(int i=0; i<Ly; i++)
		{
			// cout<<i<<","<<A[i]<<" "<<i + (Lx-1)*Ly<<","<<B[i + (Lx-1)*Ly]<<endl;
			if(A[i] != 0 && B[i + (Lx-1)*Ly] != 0 && A[i]!= B[i + (Lx-1)*Ly])
			{
				temp = B[i + (Lx-1)*Ly];
				temp1 = min(A[i], B[i + (Lx-1)*Ly]);
				// cout<<A[i]<<","<<B[i + (Lx-1)*Ly]<<endl;
				replace(B.begin(), B.end(), temp, temp1); // replaces in-place
				// if(std::find(B.begin(), B.end(), temp) != B.end()) 
				// {
				// 	cout<<"hi"<<endl;
				// }
			}
		}
	}

	if(j==1)
	{
		// cout<<"j: "<<j<<endl;
		for(int i=0; i<Ly; i++)
		{
			// cout<<i + (Lx-1)*Ly<<","<<A[i + (Lx-1)*Ly]<<" "<<i<<","<<B[i]<<endl;
			if(A[i + (Lx-1)*Ly] != 0 && B[i] != 0 && A[i + (Lx-1)*Ly]!= B[i])
			{
				temp = B[i];
				temp1 = min(A[i + (Lx-1)*Ly], B[i]);
				// cout<<A[i + (Lx-1)*Ly]<<","<<B[i]<<endl;
				replace(B.begin(), B.end(), temp, temp1); // replaces in-place
			}
		}
	}

	if(j==2)
	{
		// cout<<"j: "<<j<<endl;
		for(int i=0; i<Lx; i++)
		{
			// cout<<Lx*i<<","<<A[Lx*i]<<" "<<(Lx*i) + (Ly-1)<<","<<B[(Lx*i) + (Ly-1)]<<endl;
			if(A[Lx*i] != 0 && B[(Lx*i) + (Ly-1)] != 0 && A[Lx*i]!= B[(Lx*i) + (Ly-1)])
			{
				// cout<<A[Lx*i]<<","<<B[(Lx*i) + (Ly-1)]<<endl;
				temp = B[(Lx*i) + (Ly-1)];
				temp1 = min(A[Lx*i], B[(Lx*i) + (Ly-1)]);
				replace(B.begin(), B.end(), temp, temp1); // replaces in-place
			}
		}
	}

	if(j==3)
	{
		// cout<<"j: "<<j<<endl;
		for(int i=0; i<Lx; i++)
		{
			// cout<<(Lx*i) + (Ly-1)<<","<<A[(Lx*i) + (Ly-1)]<<" "<<Lx*i<<","<<B[Lx*i]<<endl;
			if(A[(Lx*i) + (Ly-1)] != 0 && B[Lx*i] && A[(Lx*i) + (Ly-1)]!= B[Lx*i])
			{
				// cout<<A[(Lx*i) + (Ly-1)]<<","<<B[Lx*i]<<endl;
				temp = B[Lx*i];
				temp1 = min(A[(Lx*i) + (Ly-1)], B[Lx*i]);
				replace(B.begin(), B.end(), temp, temp1); // replaces in-place
			}
		}
	}

} 

vector<int> HKPBC:: serial(vector<vector<int>> &niegh, vector<vector<int>> &total, int max_lab)
{
	for (int i=0; i<niegh.size(); i++) //niegh.size()
	{
		for (int j=0; j<4; j++)
		{
			// cout<<i<<","<<niegh[i][j]<<endl;
			Hks_pbc(total[i], total[niegh[i][j]], j);
		}

	}

	int clusters[max_lab + 1] = {0};
	int clus[max_lab];
	for(int j=0; j<total.size(); j++)
	{
		for(int i=0;i<(Lx*Ly);i++)
		{  
			clusters[total[j][i]]++;
		}
	}

	// cout<<"the orginal:"<<endl;
	// for(int j=0; j<max_lab + 1; j++)
	// {
	// 	cout<<clusters[j]<<",";
	// }

	// cout<<endl;

	copy(clusters + 1, clusters + max_lab +1, clus);

	// cout<<"cluster size"<<endl;
	// for (int i = 0; i < max_lab ; i++)
	// {
	// cout<<clus[i]<<",";
	// }
	// cout<<endl;
	
	vector<int> cluster_index(max_lab);
	std::iota(cluster_index.begin(),cluster_index.end(),0); //Initializing
	sort( cluster_index.begin(),cluster_index.end(), [&](int i,int j){return clus[i]<clus[j];} );
	std::for_each(cluster_index.begin(), cluster_index.end(), [](int& d) { d+=1;});

	// cout<<"cluster_index"<<endl;
	// for(int j=0; j<max_lab; j++)
	// {
	// 	cout<<cluster_index[j]<<",";
	// }

	// cout<<endl;

	return cluster_index;
}




