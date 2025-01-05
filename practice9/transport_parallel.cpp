#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int Nx = 20;
    const int Ny = 20;
    const double cx = 1.0;
    const double cy = 1.0;
    const double Lx = 1.0;
    const double Ly = 1.0;
    const double T = 1.0;
    const int Nt = ((Nx * Ny) / (cx + cy)); // Stability condition
    const double dx = Lx / Nx;
    const double dy = Ly / Ny;
    const double dt = T / Nt;
    double start, end;

    int k = Nx / size;
    // Same logic as previous lab, now 3d matrix
    vector<vector<vector<double>>> u_local(Nt, vector<vector<double>>(k + 1, vector<double>(Ny, 0.0)));

    if(rank==0){start = MPI_Wtime();}

    // Each processor applies initial conditions to u_local, it's better than scattering global vector
    for(int i = 0; i <k; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = (i + rank * k) * dx;
            double y = j * dy;
            u_local[0][i + 1][j] = (x > 0 && y > 0) ? 1.0 : 0.0;
        }
    }

    for(int n = 0; n < Nt - 1; ++n) {

        if(rank > 0) {
            // MPI_Send(&u_local[n][1][0], Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_local[n][0][0], Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (rank < size-1) {
            MPI_Send(&u_local[n][k][0], Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            // MPI_Recv(&u_local[n][k + 1][0], Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i <= k; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double ux = (u_local[n][i][j] - u_local[n][i - 1][j]) / dx;
                double uy = (u_local[n][i][j] - u_local[n][i][j - 1]) / dy;
                u_local[n + 1][i][j] = u_local[n][i][j] - dt * (cx * ux + cy * uy);
            }
        }
    }

    vector<vector<vector<double>>> u_global(Nt, vector<vector<double>>(Nx, vector<double>(Ny, 0.0)));

    for (int n = 0; n < Nt; ++n) {
        vector<double> buffer(Ny, 0.0);
        if (rank == 0) {
            // Store rank 0's portion
            for (int i = 1; i <= k; ++i) {
                u_global[n][i - 1] = u_local[n][i];
            }

            // Store data from other ranks
            for (int p = 1; p < size; ++p) {
                for (int i = 0; i < k; ++i) {
                    MPI_Recv(buffer.data(), Ny, MPI_DOUBLE, p, n * size + p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    u_global[n][p * k + i] = buffer;
                }
            }
        } else {
            // Other ranks send data row by row
            for (int i = 1; i <= k; ++i) {
                MPI_Send(u_local[n][i].data(), Ny, MPI_DOUBLE, 0, n * size + rank, MPI_COMM_WORLD);
            }
        }
    }

    if (rank == 0) {
        end = MPI_Wtime();
        ofstream outfile("transport_solution_2D_parallel.txt");
        if (!outfile.is_open()) {
            cerr << "Error opening file" << endl;
            MPI_Finalize();
            return 1;
        }

        for (int n = 0; n < Nt; ++n) {
            double t = n * dt;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    outfile << t << " " << i * dx << " " << j * dy << " " << u_global[n][i][j] << endl;
                }
                outfile << endl;
            }
            outfile << endl;
        }
        outfile.close();
        cout << "Saved to transport_solution_2D.txt" << endl;
        cout<<"Time taken: "<<end - start<<endl;
    }

    MPI_Finalize();
    return 0;
}
