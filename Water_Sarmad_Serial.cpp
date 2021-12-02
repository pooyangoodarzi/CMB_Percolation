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
#include <cstdlib>
#include <iostream>
#include <thread>
#include <map>
#include <string>
#include <valarray>
#include <future>
#include <iterator>

using namespace std;

#define N 50320512
#define Lx 5016
#define Ly 10032
#define Nsplit 88


vector<vector<int>> split (const std::vector<int>& v);
double cosine(vector<double> plat);
vector<int> Divisors(int n);
vector<vector<double>> water_clusters(int ii, int const& num_iter,vector<int> indexes, double total_lat, vector<double> plat, int PBS_Pos);

int main()
{
    double total_lat;
    vector<double> final_total_landmass;
    vector<double> final_big_cluster;

    // //opening data with cnpy
    // cnpy::NpyArray arr = cnpy::npy_load("/home/complex/c++/Rewrite_c++/HKPBC/pmap_22.npy");
    // double* pmap = arr.data<double>();

    // cnpy::NpyArray arr1 = cnpy::npy_load("/home/complex/c++/Rewrite_c++/HKPBC/plat_22.npy");
    // double* plat = arr1.data<double>();
    std::vector<double> pmap, plat;
    std::string line;

    std::ifstream myFile("/share/users/m_movahed/cmb_cluster/s0_c++/pmap.txt");
    while(getline(myFile, line))
    {
        std::istringstream lineStream(line);
        double first;
        lineStream >> first;
        pmap.push_back(first);
    }

    std::ifstream myFile1("/share/users/m_movahed/cmb_cluster/s0_c++/plat.txt");
    while(getline(myFile1, line))
    {
        std::istringstream lineStream(line);
        double first;
        lineStream >> first;
        plat.push_back(first);
    }
    //calculating cosine of the all the elements in plat
    total_lat = cosine(plat);
    cout<<total_lat<<endl;

    //sorting the indexes
    vector<int> indexes(N);
    std::iota(indexes.begin(),indexes.end(),0); //Initializing
    sort(indexes.begin(), indexes.end(), [&](int i,int j){return pmap[i]<pmap[j];} );

    //Defining PBS Variables
    auto PBS_index = split(indexes);
    int PBS_counter = two;
    int PBS_Pos = PBS_counter*PBS_index[0].size();
    //cout<<PBS_index[0].size()<<"_"<<PBS_index[1].size()<<endl;



    // Finding number of processes
    vector<int> divisions;
    divisions = Divisors(PBS_index[0].size());
    int num_iter = divisions[ceil(divisions.size()/2)];
    int num_proc = int(PBS_index[0].size()/num_iter);
    // int num_iter = 8;
    // int num_proc = 4;
    printf("num_proc = %d \n", num_proc);
    printf("num_iter = %d \n", num_iter);


    //using future for paralizing
    std::vector<std::future<vector<vector<double>>>> futures;

    for (int j = 0; j< num_proc; j++)
    {
        futures.push_back(std::async(std::launch::async, water_clusters, j, num_iter, indexes, total_lat, plat, PBS_Pos));
    }

    // collecting the resualt of the future
    for (int j = 0; j< num_proc; j++)
    {
        auto results = futures[j].get();
        std::copy(results[0].begin(), results[0].end(),  back_inserter(final_total_landmass));
        std::copy(results[1].begin(), results[1].end(),  back_inserter(final_big_cluster));
        printf("%d \n", j);
    }
    // std::future<vector<double>> futuer;
    // std::vector<vector<vector<double>>> futures;

    // for (int j = 0; j< 1; j++)//num_proc
    // {
    //     std::future<vector<vector<double>>> futu = std::async(std::launch::async, water_clusters, j, num_iter, indexes, total_lat, plat);
    //     auto results = futu.get();
    //     std::copy(results[0].begin(), results[0].end(),  back_inserter(final_total_landmass));
    //     std::copy(results[1].begin(), results[1].end(),  back_inserter(final_big_cluster));
    //     printf("%d \n", j);
    // }

    //saving the data in npy format
    string path1 = "/share/users/m_movahed/cmb_cluster/s0_c++/";
    string path = "/share/users/m_movahed/cmb_cluster/s0_c++/";
    string filename = "total_landmass_" + to_string(one) + "_" + to_string(PBS_counter) + ".txt";
    string filename1 = "big_cluster_" + to_string(one) + "_" + to_string(PBS_counter) + ".txt";

    // cnpy::npy_save(path + filename , &final_total_landmass[0],{N},"w");
    // cnpy::npy_save(path1 + filename1 , &final_big_cluster[0],{N},"w");

    std::ofstream fout(path + filename);
    fout.precision(17);
    std::copy(final_total_landmass.begin(), final_total_landmass.end(),std::ostream_iterator<double>(fout, "\n"));

    std::ofstream fout1(path1 + filename1);
    fout1.precision(17);
    std::copy(final_big_cluster.begin(), final_big_cluster.end(),std::ostream_iterator<double>(fout1, "\n"));

        return 0;
}

//Split Function
std::vector<std::vector<int>> split (const std::vector<int>& v)
{
    int n = v.size();
    int size_max = n / Nsplit + (n % Nsplit != 0);
    std::vector<std::vector<int>> split;
    for (int ibegin = 0; ibegin < n; ibegin += size_max) {
        int iend = ibegin + size_max;
        if (iend > n) iend = n;
        split.emplace_back (std::vector<int>(v.begin() + ibegin, v.begin() + iend));
    }
    return split;
}


//divisors function to caculate the number of processes and iterations
vector<int> Divisors(int n)
{
    // Vector to store half of the divisors
    vector<int> v;
    for (int i = 1; i <= sqrt(n); i++) {
        if (n % i == 0) {
            // check if divisors are equal
            if (n / i != i)
            {
                v.push_back(n / i);
            }
        }
    }
    return v;
}

//function of calculating all cosines of plat
double cosine(vector<double> plat)
{
    double sum = 0;
    for (int i = 0; i <N; i++)
    {
        sum = sum + abs(cos(plat[i]));
    }
    return sum;
}

//main water function
vector<vector<double>> water_clusters(int ii, int const& num_iter, vector<int> indexes, double total_lat, vector<double> plat, int PBS_Pos)
{
    int l = ii*num_iter;

    vector<double> total_landmass;
    vector<double> big_cluster;
    vector<vector<double>> result;
    double pp;
    vector<int> s;
    HKPBC object(Lx,Ly,N);


    for (int ii=l; ii< l + num_iter; ii++)
    {
        vector<int> myBoolArray(N);
        for (int m=0; m< (ii+1)+PBS_Pos; m++)
        {
            myBoolArray[indexes[m]] = 1;
        }
        s = object.HK(myBoolArray);
        vector<double> allclusters;

        if (s.size() >= 3)
        {
            for(int j=0; j < 3; j++)
            {
                valarray<int> mask ( &(myBoolArray[0]), N);
                mask[mask < s.rbegin()[j]] = 0;
                mask[mask > s.rbegin()[j]] = 0;
                mask[mask == s.rbegin()[j]] = 1;

                vector<double> masklat;
                double sum =0;
                for(int k=0; k<N; k++)
                {
                    masklat.push_back(plat[k] * mask[k]);
                    if(masklat[k] != 0)
                    {
                        sum = sum + abs(cos(masklat[k]));
                    }
                }
                allclusters.push_back(sum);
            }
        }

        if (s.size()==2)
        {
            for(int j=0; j < 2; j++)
            {
                valarray<int> mask ( &(myBoolArray[0]), N);
                mask[mask < s.rbegin()[j]] = 0;
                mask[mask > s.rbegin()[j]] = 0;
                mask[mask == s.rbegin()[j]] = 1;

                vector<double> masklat;
                double sum =0;
                for(int k=0; k<N; k++)
                {
                    masklat.push_back(plat[k] * mask[k]);
                    if(masklat[k] != 0)
                    {
                        sum = sum + abs(cos(masklat[k]));
                    }
                }

                allclusters.push_back(sum);
            }

        }
        if (s.size() == 1)
        {
            for(int j=0; j < 1; j++)
            {
                valarray<int> mask ( &(myBoolArray[0]), N);
                mask[mask < s.rbegin()[j]] = 0;
                mask[mask > s.rbegin()[j]] = 0;
                mask[mask == s.rbegin()[j]] = 1;

                vector<double> masklat;
                double sum =0;
                for(int k=0; k<N; k++)
                {
                    masklat.push_back(plat[k] * mask[k]);
                    if(masklat[k] != 0)
                    {
                        sum = sum + abs(cos(masklat[k]));
                    }
                }

                allclusters.push_back(sum);
            }
        }

        pp = ((ii+1)+PBS_Pos) / double(N);
        big_cluster.push_back(*max_element(allclusters.begin(), allclusters.end()) / total_lat);
        // cout<<"big_cluster = "<<*max_element(allclusters.begin(), allclusters.end()) / total_lat<<endl;
        total_landmass.push_back(pp);
    }

    result.push_back(total_landmass);
    result.push_back(big_cluster);

    return result;
}



