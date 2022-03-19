//this code is a copy of the hoshen kopelman code written by https://gist.github.com/tobin/909424 i change it a little bit to add the function for periodic boundary condition from left to right side of the array
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
	
void HKPBC:: uf_pbc(std::vector<int> &matrix) //this function add the periodic boundary condition
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
	  labels = static_cast<int*>(calloc(sizeof(int), n_labels));
	  labels[0] = 0;
	}
	
void HKPBC:: uf_done(void) 
	{
	  n_labels = 0;
	  free(labels);
	  labels = 0;
	}
	
vector<int> HKPBC::HK(std::vector<int> &matrix)
	{
	
		uf_initialize(Ly * Lx / 2);
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
	  

	  uf_pbc(matrix);   //applying the periodic boundary condition
	
	  int *new_labels = static_cast<int*>(calloc(sizeof(int), n_labels)); // allocate array, initialized to zero
	  
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
	  int clus;
	  for(int i=0;i<size;i++)
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

	  return cluster_index;
	}

