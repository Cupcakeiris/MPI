#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Status send_status, recv_status;
    MPI_Request send_request, recv_request;
    MPI_Comm comm = MPI_COMM_WORLD;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double start, end;
    int left = (rank - 1 + size) % size;  
    int right = (rank + 1) % size;        

    int local_sum = rank;

    if(rank == 0){
        double start = MPI_Wtime();

        MPI_Isend(&local_sum, 1, MPI_INT, right, 0, comm, &send_request);
        MPI_Wait(&send_request, &send_status);
        MPI_Recv(&local_sum, 1, MPI_INT, left, 0, comm, &recv_status);

        double end = MPI_Wtime();

        cout<<"Result: "<<local_sum<<"\tTime taken: "<<end - start;
    }
    else{
        MPI_Recv(&local_sum, 1, MPI_INT, left, 0, comm, &recv_status);
        local_sum += rank;
        MPI_Isend(&local_sum, 1, MPI_INT, right, 0, comm, &send_request);
        MPI_Wait(&send_request, &send_status);
    }

    MPI_Finalize();
    return 0;
}
