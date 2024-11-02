#include <mpi.h>
#include <iostream>
#include <vector>
using namespace std;

void print_matrix(const vector<vector<int>>& A);
int product(const vector<vector<int>>& A, const vector<int>& B) ;
void print_vector(const vector<int>& vec);

int main(int argc, char * argv[]){
    int rank, size, error_code;
    int number_of_working_ranks = 3;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Wtime(); // Time starts when MPI is initialized
    MPI_Status status;

    ///////////////////////////////////////////////////// Task 3
    int n = 7;
    int m = 6;
    vector<vector<int>> A(n, vector<int>(m)); // Ax = b
    vector<int> x(m);
    vector<int> b(n);

    if(rank==0){
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = (rand() % 10) + 1;
            }
        }
        for (int i = 0; i < m; i++) {
            x[i] = (rand() % 10) + 1;
        }
                
        cout<<"-----------Generated values:---------------\n\n";
        print_matrix(A);
        cout<<"\n";
        print_vector(x);

        // sends each row
        for (int worker_rank = 1; worker_rank < size; worker_rank++) {
            MPI_Send(A[worker_rank - 1].data(), A[0].size(), MPI_INT, worker_rank, 0, MPI_COMM_WORLD);  
            MPI_Send(x.data(), A[0].size(), MPI_INT, worker_rank, 1, MPI_COMM_WORLD); 
        }

        double start_time = MPI_Wtime();
        for (int worker_rank = 1; worker_rank < size; worker_rank++) {
            MPI_Recv(&b[worker_rank - 1], 1, MPI_FLOAT, worker_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        cout<<"------------------Product result--------------\n\n";
        double end_time = MPI_Wtime();
        print_vector(b);
        cout<<"Time spent: "<<end_time - start_time<<"\n";      
        cout<<"--------------------Analytical result-----------\n\n";
        double start_analytical = MPI_Wtime();
        product(A, x);
        double end_analytical = MPI_Wtime();
        cout<<"Time spent: "<<end_analytical - start_analytical<<"\n";
    }

    else{
        vector<int> A_row(A[0].size()); // in 1 row we have # of columns elements
        vector<int> x_copy(x.size());

        MPI_Recv(A_row.data(), A[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(x_copy.data(), A[0].size(), MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int result = 0;
        for (int i = 0; i < A[0].size(); i++) {
            result += A_row[i] * x_copy[i];
        }
        MPI_Send(&result, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }

}

void print_matrix(const vector<vector<int>>& A) {
    for (const auto& row : A) {
            for (int value : row) {
                cout << value << "\t";
            }
            cout << "\n";
        }
}

void print_vector(const vector<int>& vec){
    for (int value : vec) {
        cout << value << "\t";
    }
    cout << endl;
}

int product(const vector<vector<int>>& A, const vector<int>& B) {
    int rows_A = A.size(); 
    int cols_A = A[0].size();   
    int size_B = B.size();     

    vector<int> result(rows_A, 0);

    for (int i = 0; i < rows_A; i++) {
        result[i] = 0;
        for (int k = 0; k < cols_A; k++) {
            result[i] += A[i][k] * B[k];
        }
        cout << result[i] << "\t";
    }
    cout<<"\n";
    return 0;
}