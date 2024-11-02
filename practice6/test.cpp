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

    if(rank%2==0){
    cout<<rank<<"hi\n";}
    else{
        cout<<rank<<"bye\n";
    }

    
    MPI_Finalize();
}