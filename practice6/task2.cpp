#include <mpi.h>
#include <iostream>
#include <vector>

using namespace std;

void print_matrix(const vector<vector<double>> &A){
    for(const auto& row:A){
        for(double value:row){
            cout<<value<<"\t";
        }
        cout<<"\n";
    }
}

int main(int argc, char *argv[]){
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int n = 100; // rows in matrix A
    int m = 100; // columns in matrix A and rows in matrix B
    int p = 100; // Number of columns in matrix B

    // n has to be divisible by size
    if (n%size !=0){
        if(rank==0){
            cout<<"Number of rows n must be divisibile by number of proccesses (size)\n";
            MPI_Abort( comm, 1);
        }
    }

    int local_rows_size = n / size;

    // Initialize matrices
    vector<vector<double>> A(n, vector<double>(m));
    vector<vector<double>> B(m, vector<double>(p));
    vector<vector<double>> result(n, vector<double>(p, 0));

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = (rand() % 10) + 1;
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                B[i][j] = (rand() % 10) + 1;
            }
        }

        cout << "-----------Generated Matrix A---------------\n\n";
        print_matrix(A);
        cout << "\n-----------Generated Matrix B---------------\n\n";
        print_matrix(B);
        cout << "\n";
    }

    double start = MPI_Wtime();
    for(int i = 0; i < m; i++){ // Bcasting B's rows
        MPI_Bcast(B[i].data(), p, MPI_DOUBLE, 0, comm);
    }

    // flatten A to scatter
    vector<double> A_flat(n * m);
    if(rank == 0){
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                A_flat[i*m + j] = A[i][j];
            }
        }
    }

    vector<double> local_row(local_rows_size * m);
    MPI_Scatter( A_flat.data(), local_rows_size * m, MPI_DOUBLE, 
                local_row.data(), local_rows_size * m, MPI_DOUBLE, 0, comm);

    // local matrix multiplication
    vector<vector<double>> local_matrix(local_rows_size, vector<double>(p, 0));
    for (int i = 0; i<local_rows_size; i++){
        for(int j=0; j<p; j++){
            for(int k=0; k<m; k++){
                local_matrix[i][j] += local_row[i*m + k] * B[k][j];
        }
        }
    }

    // flatten local result
    vector<double> flat_local_matrix(local_rows_size * p);
    for(int i = 0; i<local_rows_size; i++){
        for(int j = 0; j<p; j++){
            flat_local_matrix[i * p + j] = local_matrix[i][j];
        }
    }

    // gather results
    vector<double> flat_result(n*p);
    MPI_Gather( flat_local_matrix.data(), local_rows_size*p , MPI_DOUBLE,
                flat_result.data(), local_rows_size*p, MPI_DOUBLE , 0, comm);

    // unflatten results
    if(rank == 0){
        for(int i =0; i<n; i++){
            for(int j=0; j<p; j++){
                result[i][j] = flat_result[i*p+j];
            }
        }

        double end = MPI_Wtime();

        cout << "------------------Result-----------------\n\n";
        print_matrix(result);
        cout<<"\nTime taken: "<<end - start<<"\n";
    }

    MPI_Finalize();
}