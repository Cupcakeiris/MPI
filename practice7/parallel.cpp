#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <mpi.h>
using namespace std;

void save_to_file(const vector<float> &arr) {
    ofstream outfile("newP_parallelized.txt");
    for (size_t i = 0; i < arr.size(); i++) {
        outfile << arr[i] << "\n";
    }
    outfile.close();
}

float f(float x) {
    return 2 * exp(x) * cos(x);
}

int main(int argc, char *argv[]) {
    int Nx = 1000;
    float L = 1.0;
    float dx = L / (Nx - 1);
    float tol = 1e-6;
    float max_err = 100.0;
    float local_max_err;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double start, end;

    vector<float> oldP(Nx, 0);
    vector<float> newP(Nx, 0);

    int k = Nx / size;
    int i1 = rank * k;

    vector<float> local_oldP(k + 2, 0);
    vector<float> local_newP(k + 2, 0);

    MPI_Scatter(oldP.data(), k, MPI_FLOAT, &local_oldP[1], k, MPI_FLOAT, 0, comm);
    MPI_Scatter(newP.data(), k, MPI_FLOAT, &local_newP[1], k, MPI_FLOAT, 0, comm);

    // Boundary conditions
    if (rank == 0) {
        local_oldP[0] = 1.0;
        local_newP[0] = 1.0;
        start = MPI_Wtime();
    }
    if (rank == size - 1) {
        local_oldP[k] = 0.0;
        local_newP[k] = 0.0;
    }

    while (max_err > tol) {
        local_max_err = 0.0;

        // Send and receive ghost nodes for boundaries
        if (rank < size - 1) {
            MPI_Send(&local_oldP[k], 1, MPI_FLOAT, rank + 1, 0, comm);
            MPI_Recv(&local_oldP[k + 1], 1, MPI_FLOAT, rank + 1, 0, comm, MPI_STATUS_IGNORE);
        }
        if (rank > 0) {
            MPI_Send(&local_oldP[1], 1, MPI_FLOAT, rank - 1, 0, comm);
            MPI_Recv(&local_oldP[0], 1, MPI_FLOAT, rank - 1, 0, comm, MPI_STATUS_IGNORE);
        }

        // Perform computation excluding ghost nodes
        for (int i = 1; i < k + 1; i++) {
            float x = (i1 + i - 1) * dx;  // for global x
            local_newP[i] = 0.5 * (local_oldP[i - 1] + local_oldP[i + 1] + dx * dx * f(x));
            local_max_err = max(local_max_err, abs(local_newP[i] - local_oldP[i]));
        }

        local_oldP = local_newP;

        // Compute the global maximum error
        MPI_Allreduce(&local_max_err, &max_err, 1, MPI_FLOAT, MPI_MAX, comm);

    }
        // Gather results excluding ghost nodes
        MPI_Gather(&local_newP[1], k, MPI_FLOAT, newP.data(), k, MPI_FLOAT, 0, comm);

    if (rank == 0) {
        save_to_file(newP);
        cout << "Data saved to newP_parallelized.txt" << endl;
        end = MPI_Wtime();
        cout<<"Time taken: "<<end - start<<endl;
    }

    MPI_Finalize();
    return 0;
}
