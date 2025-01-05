#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;
// Same logic as Lax

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int Nx = 100; 
    const int Nt = 100; 
    const double L = 1.0; 
    const double T = 1.0; 
    const double dx = L / (Nx - 1); 
    const double dt = T / Nt; 
    const double c = 1.0; 
    const double alpha = c * dt / (dx * dx);
    double start, end;

    int k = Nx / size;  

    vector<vector<double>> u_local(Nt, vector<double>(k + 2, 0.0));

    vector<double> u(Nx, 0.0);
    if (rank == 0) {
        for (int i = 0; i < Nx; ++i) {
            u[i] = (i * dx > 0) ? 1.0 : 0.0; 
        }
        start = MPI_Wtime();
    }

    MPI_Scatter(u.data(), k, MPI_DOUBLE, &u_local[0][1], k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int n = 0; n < Nt - 1; ++n) {

        if (rank > 0) {
            MPI_Send(&u_local[n][1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_local[n][0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(&u_local[n][k], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_local[n][k + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i <= k; ++i) {
            u_local[n + 1][i] = u_local[n][i] - alpha * (u_local[n][i] - u_local[n][i - 1]);
        }
    }

    vector<vector<double>> u_global(Nt, vector<double>(Nx, 0.0));
    for (int n = 0; n < Nt; ++n) {
        if (rank == 0) {
            // Store Rank 0's portion
            for (int i = 1; i <= k; ++i) {
                u_global[n][i - 1] = u_local[n][i];
            }

            // Store data from other ranks
            for (int p = 1; p < size; ++p) {
                vector<double> buffer(k, 0.0);
                MPI_Recv(buffer.data(), k, MPI_DOUBLE, p, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < k; ++i) {
                    u_global[n][p * k + i] = buffer[i];
                }
            }
        } else {
            MPI_Send(u_local[n].data(), k, MPI_DOUBLE, 0, n, MPI_COMM_WORLD);
        }
    }


    if (rank == 0) {
        end = MPI_Wtime();
        ofstream outfile("numerical_solution_back.txt");
        if (!outfile.is_open()) {
            cerr << "Error opening file" << endl;
            MPI_Finalize();
            return 1;
        }

        for (int n = 0; n < Nt; ++n) {
            double t = n * dt;
            for (int i = 0; i < Nx; ++i) {
                outfile << t << " " << i * dx << " " << u_global[n][i] << endl;
            }
            outfile << endl;
        }

        outfile.close();
        cout << "Saved to numerical_solution_back.txt" << endl;
        cout << "Time taken: " << end - start << " seconds" << endl;
    }

    MPI_Finalize();
    return 0;
}
