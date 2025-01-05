#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int Nx = 1000;  
    const int Nt = 2000;
    const double L = 10.0; 
    const double T = 6.0; 
    const double dx = L / (Nx - 1); 
    const double c = 2.0; 
    const double dt = T / Nt;
    const double stable = dx / c;
    const double alpha = c * dt / dx;
    
    if (rank == 0 && dt > stable) {
        cout << "Scheme is not stable\n";
        MPI_Abort(MPI_COMM_WORLD, 1); // Abort if unstable
        return 1;
    }

    int k = Nx / size;  

    vector<vector<double>> u_local(Nt, vector<double>(k + 2, 0.0));

    vector<double> u(Nx);
    if (rank == 0) {
        for (int i = 0; i < Nx; ++i) {
            u[i] = cos(i * dx) * cos(i * dx);
        }
    }

    // Scatter initial conditions
    MPI_Scatter(u.data() + rank * k, k, MPI_DOUBLE, &u_local[0][1], k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double start = MPI_Wtime();

    for (int n = 0; n < Nt - 1; n++) {
        if (rank == 0) {
            u_local[n][0] = 1.0;
        }

        // Communicate boundary values
        if (rank > 0) {
            MPI_Send(&u_local[n][1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_local[n][0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(&u_local[n][k], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_local[n][k + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Lax - Friedrichs method
        for (int i = 1; i <= k; ++i) {
            u_local[n + 1][i] = 0.5 * (u_local[n][i - 1] + u_local[n][i + 1]) - 
                                 alpha * 0.5 * (u_local[n][i + 1] - u_local[n][i - 1]);

        }
    }

    vector<double> u_global(Nx);
    vector<string> output;

    for (int n = 0; n < Nt; ++n) {
        if (rank == 0) {
            // Store Rank 0's portion
            for (int i = 1; i <= k; ++i) {
                u_global[i - 1] = u_local[n][i];
            }

            // Gather data from other ranks
            for (int p = 1; p < size; ++p) {
                vector<double> buffer(k);
                MPI_Recv(buffer.data(), k, MPI_DOUBLE, p, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < k; ++i) {
                    u_global[p * k + i] = buffer[i];
                }
            }

            // Store results for this time step
            double t = n * dt;
            string timestep_output;
            for (int i = 0; i < Nx; ++i) {
                timestep_output += to_string(t) + " " + to_string(i * dx) + " " + to_string(u_global[i]) + "\n";
            }
            output.push_back(timestep_output);

        } else {
            MPI_Send(u_local[n].data() + 1, k, MPI_DOUBLE, 0, n, MPI_COMM_WORLD);
        }
    }

    if (rank == 0) {
        // Write results to file
        ofstream outfile("numerical_solution_lax.txt");
        if (!outfile.is_open()) {
            cerr << "Error opening file" << endl;
            MPI_Finalize();
            return 1;
        }

        for (const auto &timestep_output : output) {
            outfile << timestep_output << endl;
        }

        outfile.close();
        double end = MPI_Wtime();
        cout << "Saved to numerical_solution_lax.txt" << endl;
        cout << "Time taken: " << end - start << " seconds" << endl;
    }

    MPI_Finalize();
    return 0;
}
