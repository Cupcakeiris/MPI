#include <mpi.h>
#include <stdio.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
    int rank, size;
    int tag1 = 1337;
    int tag2 = 3117;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char ping = 'A';
    char pong = 'a';

    int i = 0;
    while (i < 10) {
        cout<<i<<endl;

        if (rank == 0) {
            MPI_Send(&ping, 1, MPI_CHAR, 1, tag1, MPI_COMM_WORLD);
            // printf("Process 0 sent '%c' to 1\n", ping);
            cout<<"Process 0 sent "<<ping<<" to 1\n";
            MPI_Recv(&pong, 1, MPI_CHAR, 1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("Process 0 received '%c' from 1\n", pong);
            cout<<"Process 0 received "<<pong<<" from 1\n";
            ping++;
        } else if (rank == 1) {
            MPI_Recv(&ping, 1, MPI_CHAR, 0, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("Process 1 received '%c' from 0\n", ping);
            cout<<"Process 1 received "<<ping<<" from 0\n";
            MPI_Send(&pong, 1, MPI_CHAR, 0, tag2, MPI_COMM_WORLD);
            // printf("Process 1 sent '%c' to 0\n", pong);
            cout<<"Process 1 sent "<<pong<<" to 0\n";
            pong++;
        }
        i++;
    }

    MPI_Finalize();
    return 0;
}

