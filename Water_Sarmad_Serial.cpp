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
#include <cstdlib>
#include <iostream>
#include <thread>
#include <map>
#include <string>
#include <valarray>
#include <future>
#include <iterator>

using namespace std;


#define Lx 90
#define Ly 90
#define Nsplit 1
#define patch_count 8
#define N Lx*Ly*patch_count
#define PI 3.14159265358979323846


vector<vector<int>> split (const vector<int>& v, int Num);
vector<int> flatten(const vector<vector<int>> &orig);
double cosine(vector<double> plat);
vector<int> Divisors(int n);
vector<vector<double>> water_clusters(int ii, int const& num_iter,vector<int> indexes, double total_lat, vector<double> plat, int PBS_Pos, vector<vector<int>> niegh);

int main()
{
    double total_lat;
    vector<double> final_total_landmass;
    vector<double> final_big_cluster;
    // vector<vector<int>> niegh;
    vector<double> pmap, plat, pmap_tot;
    string line;
    vector<int> temp;


    for(int files = 0; files < patch_count ;files++)
    {   
        ostringstream filename;
        filename << "/home/complex/HK_Serial/Data/" << files << ".txt";
        ifstream input(filename.str());
        while(getline(input, line))
        {
            std::istringstream lineStream(line);
            double first;
            lineStream >> first;
            pmap.push_back(first);
            pmap_tot.push_back(first);

        } 
    }



    string line1;
    for(int files = 0; files < patch_count ;files++)
    {   
        ostringstream filename;
        filename << "/home/complex/HK_Serial/Lat/" << files << ".txt";
        ifstream input(filename.str());
        while(getline(input, line1))
        {
            std::istringstream lineStream(line1);
            double first;
            lineStream >> first;
            plat.push_back(first);
        } 

    }

    std::string line2;
    std::ifstream myFile4("/home/complex/HK_Serial/nieghearth.txt"); 
    while(getline(myFile4, line2))
    {
        std::istringstream lineStream(line2);
        int first1;
        lineStream >> first1;
        temp.push_back(first1);
    }
    cout<<"hi"<<endl;
    auto niegh = split(temp, patch_count);

    //calculating cosine of the all the elements in plat
    total_lat = cosine(plat);
    cout<<total_lat<<endl;

    //sorting the indexes
    vector<int> indexes(patch_count*Lx*Ly);
    std::iota(indexes.begin(),indexes.end(),0); //Initializing
    sort(indexes.begin(), indexes.end(), [&](int i,int j){return pmap_tot[i]<pmap_tot[j];} );
    // for (size_t i = 0; i < indexes.size(); i++)
    // {
    //     cout<<indexes[i]<<" ";
    // }
    // cout<<endl;
    
    //Defining PBS Variables
    auto PBS_index = split(indexes, Nsplit);
    int PBS_counter = 0;
    int PBS_Pos = PBS_counter*PBS_index[0].size();
    // cout<<PBS_index[0].size()<<endl;



    // Finding number of processes
    vector<int> divisions;
    divisions = Divisors(PBS_index[0].size());
    int num_iter = divisions[ceil(divisions.size()/2)];
    int num_proc = int(PBS_index[0].size()/num_iter);
    // int num_iter = N;
    // int num_proc = 1;
    printf("num_proc = %d \n", num_proc);
    printf("num_iter = %d \n", num_iter);


    //using future for paralizing
    std::vector<std::future<vector<vector<double>>>> futures;

    for (int j = 0; j< num_proc; j++)
    {
        futures.push_back(std::async(std::launch::async, water_clusters, j, num_iter, indexes, total_lat, plat, PBS_Pos, niegh));
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
    //     std::future<vector<vector<double>>> futu = std::async(std::launch::async, water_clusters,  j, 1, indexes, total_lat, plat, PBS_Pos, niegh);
    //     auto results = futu.get();
    //     std::copy(results[0].begin(), results[0].end(),  back_inserter(final_total_landmass));
    //     std::copy(results[1].begin(), results[1].end(),  back_inserter(final_big_cluster));
    //     printf("%d \n", j);
    // }

    //saving the data in npy format
    string path1 = "/home/complex/HK_Serial/";
    string path = "/home/complex/HK_Serial/";
    string filename = "total_landmass_" + to_string(8) + "_" + to_string(PBS_counter) + ".txt";
    string filename1 = "big_cluster_" + to_string(8) + "_" + to_string(PBS_counter) + ".txt";


    std::ofstream fout(path + filename);
    fout.precision(17);
    std::copy(final_total_landmass.begin(), final_total_landmass.end(),std::ostream_iterator<double>(fout, "\n"));

    std::ofstream fout1(path1 + filename1);
    fout1.precision(17);
    std::copy(final_big_cluster.begin(), final_big_cluster.end(),std::ostream_iterator<double>(fout1, "\n"));

    return 0;
}

//flatten Function
vector<int> flatten(const vector<vector<int>> &orig)
{   
    vector<int> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());                                                                                         
    return ret;
} 

//Split Function
vector<vector<int>> split(const std::vector<int>& v, int Num)
{
    int n = v.size();
    int size_max = int(n/Num) ;
    vector<vector<int>> split;
    for (int ibegin = 0; ibegin < n; ibegin += size_max) 
    {
        int iend = ibegin + size_max;
        if (iend > n) iend = n;
        split.emplace_back (vector<int>(v.begin() + ibegin, v.begin() + iend));
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
    for (int i = 0; i < plat.size(); i++)
    {
        sum = sum + abs(cos(plat[i]));
    }
    return sum;
}

//main water function
vector<vector<double>> water_clusters(int ii, int const& num_iter, vector<int> indexes, double total_lat, vector<double> plat, int PBS_Pos, vector<vector<int>> niegh)
{
    int l = ii*num_iter;
    vector<double> total_landmass;
    vector<double> big_cluster;
    vector<vector<double>> result;
    double pp;
    HKPBC object(Lx,Ly,Lx*Ly);
    vector<int> s;
    for (int ii=l; ii< l + num_iter; ii++) 
    {
        vector<int> myBoolArray(N);
        for (int m=0; m< (ii+1)+PBS_Pos; m++) 
        {
            myBoolArray[indexes[m]] = 1;
        }

        auto newarr = split(myBoolArray, patch_count);

        // cout<<"hi"<<endl;
        s = object.HK(newarr, niegh);
        // cout<<"hihi"<<endl;
        auto myvector = flatten(newarr);

        // cout<<"size: "<<myvector.size()<<endl;
        // for (int i = 0; i < s.size(); i++)
        // {
        //     cout<<s[i]<<",";
        // }
        // cout<<endl;
        

        vector<double> allclusters;

        if (s.size() >= 3)
        {
            for(int j=0; j < 3; j++)
            {
                valarray<int> mask ( &(myvector[0]), N);
                mask[mask < s.rbegin()[j]] = 0;
                mask[mask > s.rbegin()[j]] = 0;
                mask[mask == s.rbegin()[j]] = 1;

                vector<double> masklat;
                for(int k=0; k<N; k++)
                {
                    if(mask[k] != 0)
                    {
                        masklat.push_back(plat[k] * mask[k]);
                    }
                }
                allclusters.push_back(cosine(masklat));
            }
        }

        if (s.size()==2)
        {
            for(int j=0; j < 2; j++)
            {
                valarray<int> mask ( &(myvector[0]), N);
                mask[mask < s.rbegin()[j]] = 0;
                mask[mask > s.rbegin()[j]] = 0;
                mask[mask == s.rbegin()[j]] = 1;

                vector<double> masklat;
                for(int k=0; k<N; k++)
                {
                    if(mask[k] != 0)
                    {
                        masklat.push_back(plat[k] * mask[k]);
                    }
                }
                allclusters.push_back(cosine(masklat));
            }

        }
        if (s.size() == 1)
        {
            for(int j=0; j < 1; j++)
            {
                valarray<int> mask ( &(myvector[0]), N);
                mask[mask < s.rbegin()[j]] = 0;
                mask[mask > s.rbegin()[j]] = 0;
                mask[mask == s.rbegin()[j]] = 1;

                vector<double> masklat;
                double sum =0;
                for(int k=0; k<N; k++)
                {
                    if(mask[k] != 0)
                    {
                        masklat.push_back(plat[k] * mask[k]);
                    }
                }
                allclusters.push_back(cosine(masklat));
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



