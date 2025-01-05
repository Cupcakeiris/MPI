#include <iostream>
#include <vector>
#include <mpi.h>
#include <bits/stdc++.h> 
using namespace std;

void randomizer(vector<double>& rand_arr){
    srand(time(0)); 

    for (auto& val : rand_arr) {
        val = rand() % 1000;
    }
}

void print(vector<double>& arr){
        for(auto val : arr){
        cout<<val<<" ";
    }
    cout<<"\n";
}

void min_diff(void* in_, void* inout_, int* count, MPI_Datatype* data_type){
    double* in = (double*)in_;
    double* inout = (double*)inout_;
    
    for(int i = 0; i < *count; i++){
        inout[i] = min(in[i], inout[i]);
        cout<<in[i]<<"\t"<<inout[i]<<endl;
    }
}

int main(int argc, char * argv[]){
    MPI_Init(&argc , &argv);
    double start, end;
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    int root = 0;

    int k = 100;
    vector<double> rand_arr(size * k, 0.0);
    vector<double> local_arr(k, 0.0);

    if(rank == root){
        start = MPI_Wtime();
        randomizer(rand_arr);
    }

    MPI_Scatter(rand_arr.data(), k, MPI_DOUBLE, local_arr.data(), k, MPI_DOUBLE, root, comm);
    
    // cout<<"R: "<<rank<<"\t";
    // print(local_arr);

    double local_max = *max_element(local_arr.begin(), local_arr.end());
    double local_min = *min_element(local_arr.begin(), local_arr.end());
    double local_diff = local_max - local_min;

    cout << "Rank " << rank << " local max: " << local_max
         << ", local min: " << local_min
         << ", local diff: " << local_diff << endl;

    MPI_Op handler;
    MPI_Op_create(&min_diff, 1, &handler);

    double global_min;
    MPI_Reduce(&local_diff, &global_min, 1, MPI_DOUBLE, handler, root, comm);

    if(rank == 0){
        end = MPI_Wtime();
        cout<<"Global min difference: "<<global_min<<endl;
        cout<<"Time taken: "<<end - start<<endl;

        ofstream outfile("execution_times.txt", ios::app);
        outfile << size << "\t" << end - start << "\n";
        outfile.close();
    }


    MPI_Finalize();
}