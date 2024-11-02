#include <iostream>
#include <mpi.h>
using namespace std;

int main(int argc, char * argv[]){
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Wtime();
    MPI_Request send_request, recv_request;
    MPI_Status status;
    double start, end;

    int partial_sum = rank;
    int recv_sum = 0;
    int left = (rank - 1 + size) % size; 
    int right = (rank + 1) % size;        

//------------------------------------ Non-blocking
    start = MPI_Wtime();
    MPI_Isend(&partial_sum, 1, MPI_INT, right, 0, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&recv_sum, 1, MPI_INT, left, 0, MPI_COMM_WORLD, &recv_request);
    
    MPI_Wait(&send_request, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_request, &status);
    partial_sum += recv_sum;

    end = MPI_Wtime();
    if (rank == 0) {
        cout << "Isend\tResult: " << partial_sum << "\t Time taken: " << end - start <<"\n";
    }

//------------------------------------ Ssend & Irecv
    partial_sum = rank;
    recv_sum = 0;
    start = MPI_Wtime();
    MPI_Irecv(&recv_sum, 1, MPI_INT, left, 0, MPI_COMM_WORLD, &recv_request);
    MPI_Ssend(&partial_sum, 1, MPI_INT, right, 0, MPI_COMM_WORLD); // somehow changing Ssend and Irecv's orders affect so much to the code

    MPI_Wait(&recv_request, &status);
    partial_sum += recv_sum;

    end = MPI_Wtime();
    if (rank == 0) {
        cout << "Ssend\tResult: " << partial_sum << "\t Time taken: " << end - start <<"\n";
    }

//------------------------------------ Issend & Irecv
    partial_sum = rank;
    recv_sum = 0;
    start = MPI_Wtime();
    partial_sum += recv_sum;
    MPI_Issend(&partial_sum, 1, MPI_INT, right, 0, MPI_COMM_WORLD, &send_request);
    
    MPI_Irecv(&recv_sum, 1, MPI_INT, left, 0, MPI_COMM_WORLD, &recv_request);

    MPI_Wait(&recv_request, &status);
    

    end = MPI_Wtime();
    if (rank == 0) {
        cout << "Issend\tResult: " << partial_sum << "\t Time taken: " << end - start <<"\n";
        
    }


    MPI_Finalize();
    return 0;
}
