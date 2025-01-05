#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <mpi.h>

using namespace std;

const double Lx = 1.0, Ly = 1.0;
const int Nx = 50, Ny = 50;
const double dx = Lx / Nx, dy = Ly / Ny;
const double t = 1.0;
const int Nt = 2000;
const double dt = t / Nt;
const double Re = 10.0;

// Initialize conditions
double initial_u(double x, double y) { return 0.0; }
double initial_v(double x, double y) { return 0.0; }

double norm(double x, double y) {
    return sqrt(x * x + y * y) + 1e-10;
}


void save_results(vector<vector<double>> &U, vector<vector<double>> &V, 
                  int local_Nx, int local_Ny, int global_x_start, int global_y_start, 
                  int rank, int size, int dims[], MPI_Comm cart_comm) {
    vector<double> local_U_flat(local_Nx * local_Ny);
    vector<double> local_V_flat(local_Nx * local_Ny);
    vector<double> global_U, global_V;

    for (int i = 1; i <= local_Nx; ++i) {
        for (int j = 1; j <= local_Ny; ++j) {
            int idx = (i - 1) * local_Ny + (j - 1);
            local_U_flat[idx] = U[i][j];
            local_V_flat[idx] = V[i][j];
        }
    }

    if (rank == 0) {
        global_U.resize((Nx + 1) * (Ny + 1));
        global_V.resize((Nx + 1) * (Ny + 1));
    }

    MPI_Gather(local_U_flat.data(), local_Nx * local_Ny, MPI_DOUBLE,
               global_U.data(), local_Nx * local_Ny, MPI_DOUBLE, 
               0, cart_comm);

    MPI_Gather(local_V_flat.data(), local_Nx * local_Ny, MPI_DOUBLE,
               global_V.data(), local_Nx * local_Ny, MPI_DOUBLE, 
               0, cart_comm);

    if (rank == 0) {
        ofstream file("parallel_burgers_solution.txt");

        vector<vector<double>> U_full(Nx + 1, vector<double>(Ny + 1, 0.0));
        vector<vector<double>> V_full(Nx + 1, vector<double>(Ny + 1, 0.0));

        for (int p = 0; p < size; ++p) {
            int coords[2];
            MPI_Cart_coords(cart_comm, p, 2, coords);

            int x_offset = coords[0] * (Nx / dims[0]);
            int y_offset = coords[1] * (Ny / dims[1]);

            for (int i = 0; i < local_Nx; ++i) {
                for (int j = 0; j < local_Ny; ++j) {
                    int idx = i * local_Ny + j;
                    U_full[x_offset + i][y_offset + j] = global_U[p * local_Nx * local_Ny + idx];
                    V_full[x_offset + i][y_offset + j] = global_V[p * local_Nx * local_Ny + idx];
                }
            }
        }

        // Write results to file
        for (int i = 0; i <= Nx; ++i) {
            for (int j = 0; j <= Ny; ++j) {
                double x = i * dx, y = j * dy;
                file << x << " " << y << " " << U_full[i][j] << " " << V_full[i][j] << "\n";
            }
            file << "\n";
        }
        file.close();
        cout << "Results saved to 'parallel_burgers_solution.txt'" << endl;
    }
}

// Exchange ghost layers
void exchange_boundaries(vector<vector<double>> &U, vector<vector<double>> &V,
                         int rank, int up, int down, int left, int right,
                         int local_Nx, int local_Ny, MPI_Comm cart_comm) {
    MPI_Status status;

    vector<double> send_up(local_Ny), recv_down(local_Ny);
    vector<double> send_down(local_Ny), recv_up(local_Ny);
    vector<double> send_left(local_Nx), recv_right(local_Nx);
    vector<double> send_right(local_Nx), recv_left(local_Nx);

    // --- Exchange U matrix ---
    // Send/receive rows (up/down neighbors)
    if (up != MPI_PROC_NULL) {
        for (int j = 1; j <= local_Ny; ++j) send_up[j - 1] = U[1][j];
        MPI_Sendrecv(send_up.data(), local_Ny, MPI_DOUBLE, up, 0,
                     recv_down.data(), local_Ny, MPI_DOUBLE, down, 0,
                     cart_comm, &status);
        for (int j = 1; j <= local_Ny; ++j) U[local_Nx + 1][j] = recv_down[j - 1];
    }

    if (down != MPI_PROC_NULL) {
        for (int j = 1; j <= local_Ny; ++j) send_down[j - 1] = U[local_Nx][j];
        MPI_Sendrecv(send_down.data(), local_Ny, MPI_DOUBLE, down, 1,
                     recv_up.data(), local_Ny, MPI_DOUBLE, up, 1,
                     cart_comm, &status);
        for (int j = 1; j <= local_Ny; ++j) U[0][j] = recv_up[j - 1];
    }

    // Send/receive columns (left/right neighbors)
    if (left != MPI_PROC_NULL) {
        for (int i = 1; i <= local_Nx; ++i) send_left[i - 1] = U[i][1];
        MPI_Sendrecv(send_left.data(), local_Nx, MPI_DOUBLE, left, 2,
                     recv_right.data(), local_Nx, MPI_DOUBLE, right, 2,
                     cart_comm, &status);
        for (int i = 1; i <= local_Nx; ++i) U[i][local_Ny + 1] = recv_right[i - 1];
    }

    if (right != MPI_PROC_NULL) {
        for (int i = 1; i <= local_Nx; ++i) send_right[i - 1] = U[i][local_Ny];
        MPI_Sendrecv(send_right.data(), local_Nx, MPI_DOUBLE, right, 3,
                     recv_left.data(), local_Nx, MPI_DOUBLE, left, 3,
                     cart_comm, &status);
        for (int i = 1; i <= local_Nx; ++i) U[i][0] = recv_left[i - 1];
    }

    // --- Exchange V matrix ---
    if (up != MPI_PROC_NULL) {
        for (int j = 1; j <= local_Ny; ++j) send_up[j - 1] = V[1][j];
        MPI_Sendrecv(send_up.data(), local_Ny, MPI_DOUBLE, up, 4,
                     recv_down.data(), local_Ny, MPI_DOUBLE, down, 4,
                     cart_comm, &status);
        for (int j = 1; j <= local_Ny; ++j) V[local_Nx + 1][j] = recv_down[j - 1];
    }

    if (down != MPI_PROC_NULL) {
        for (int j = 1; j <= local_Ny; ++j) send_down[j - 1] = V[local_Nx][j];
        MPI_Sendrecv(send_down.data(), local_Ny, MPI_DOUBLE, down, 5,
                     recv_up.data(), local_Ny, MPI_DOUBLE, up, 5,
                     cart_comm, &status);
        for (int j = 1; j <= local_Ny; ++j) V[0][j] = recv_up[j - 1];
    }

    if (left != MPI_PROC_NULL) {
        for (int i = 1; i <= local_Nx; ++i) send_left[i - 1] = V[i][1];
        MPI_Sendrecv(send_left.data(), local_Nx, MPI_DOUBLE, left, 6,
                     recv_right.data(), local_Nx, MPI_DOUBLE, right, 6,
                     cart_comm, &status);
        for (int i = 1; i <= local_Nx; ++i) V[i][local_Ny + 1] = recv_right[i - 1];
    }

    if (right != MPI_PROC_NULL) {
        for (int i = 1; i <= local_Nx; ++i) send_right[i - 1] = V[i][local_Ny];
        MPI_Sendrecv(send_right.data(), local_Nx, MPI_DOUBLE, right, 7,
                     recv_left.data(), local_Nx, MPI_DOUBLE, left, 7,
                     cart_comm, &status);
        for (int i = 1; i <= local_Nx; ++i) V[i][0] = recv_left[i - 1];
    }
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    int periods[2] = {0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);

    int coords[2], up, down, left, right;
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    MPI_Cart_shift(cart_comm, 0, 1, &left, &right);
    MPI_Cart_shift(cart_comm, 1, 1, &down, &up);

    int local_Nx = Nx / dims[0];
    int local_Ny = Ny / dims[1];

    int global_x_start = coords[1] * local_Nx;
    int global_y_start = coords[0] * local_Ny;

    vector<vector<double>> U(local_Nx + 2, vector<double>(local_Ny + 2, 0.0));
    vector<vector<double>> V(local_Nx + 2, vector<double>(local_Ny + 2, 0.0));
    
    for (int i = 1; i <= local_Nx; ++i) {
        for (int j = 1; j <= local_Ny; ++j) {
            double x = (global_x_start + i - 1) * dx;
            double y = (global_y_start + j - 1) * dy;
            U[i][j] = initial_u(x, y);
            V[i][j] = initial_v(x, y);
        }
    }
    vector<vector<double>> U_new = U, V_new = V;


        // Apply boundary conditions
            if (coords[1] == 0) { // x = 0
                for (int j = 1; j <= local_Ny; ++j) {
                    U[0][j] = 0.0;
                    U_new[0][j] = 0.0;

                    V[0][j] = 0.0;
                    V_new[0][j] = 0.0;
                }
            }

            if (coords[1] == dims[0] - 1) { // x = 1
                for (int j = 1; j <= local_Ny; ++j) {
                    double y = (global_y_start + j - 1) * dy;
                    if (y > 0.7 && y < 1.0) {
                        U[local_Nx + 1][j] = -1.0;
                        U_new[local_Nx + 1][j] = -1.0;

                        V[local_Nx + 1][j] = 0.0;
                        V_new[local_Nx + 1][j] = 0.0;
                    } else {
                        U[local_Nx + 1][j] = 0.0;
                        U_new[local_Nx + 1][j] = 0.0;

                        V[local_Nx + 1][j] = 0.0;
                        V_new[local_Nx + 1][j] = 0.0;
                    }
                }
            }

            if (coords[0] == 0) { // y = 0
                for (int i = 1; i <= local_Nx; ++i) {
                    double x = (global_x_start + i - 1) * dx;
                    if (x < 0.3) {
                        U[i][0] = 0.0;
                        U_new[i][0] = 0.0;

                        V[i][0] = 0.0;
                        V_new[i][0] = 0.0;
                    } else if (x >= 0.3 && x <= 0.6) {
                        U[i][0] = U[i][1]; // Neumann
                        U_new[i][0] = U_new[i][1];

                        V[i][0] = V[i][1];
                        V_new[i][0] = V_new[i][1];
                    } else {
                        U[i][0] = 0.0;
                        U_new[i][0] = 0.0;
                        V[i][0] = 0.0;
                        V_new[i][0] = 0.0;
                    }
                }
            }

            if (coords[0] == dims[1] - 1) { // y = 1
                for (int i = 1; i <= local_Nx; ++i) {
                    double x = (global_x_start + i - 1) * dx;
                    if (x < 0.3) {
                        U[i][local_Ny + 1] = U[i][local_Ny]; // Neumann
                        U_new[i][local_Ny + 1] = U_new[i][local_Ny];

                        V[i][local_Ny + 1] = V[i][local_Ny];
                        V_new[i][local_Ny + 1] = V_new[i][local_Ny];
                    } else {
                        U[i][local_Ny + 1] = 0.0;
                        U_new[i][local_Ny + 1] = 0.0;

                        V[i][local_Ny + 1] = 0.0;
                        V_new[i][local_Ny + 1] = 0.0;
                    }
                }
            }
       



    // Main loop
    for (int n = 0; n < Nt; ++n) {
        exchange_boundaries(U, V, rank, up, down, left, right, local_Nx, local_Ny, cart_comm);

        // Use explicit time-stepping
        for (int i = 1; i <= local_Nx; ++i) {
            for (int j = 1; j <= local_Ny; ++j) {
                double x = (global_x_start + i - 1) * dx; // Correct x position
                double y = (global_y_start + j - 1) * dy; // Correct y position
                double norm_n = norm(x, y);

                double dUdn = (U[i + 1][j] - U[i - 1][j]) * (0.5 * x / norm_n) +
                              (U[i][j + 1] - U[i][j - 1]) * (0.5 * y / norm_n);
                double dVdn = (V[i + 1][j] - V[i - 1][j]) * (0.5 * x / norm_n) +
                              (V[i][j + 1] - V[i][j - 1]) * (0.5 * y / norm_n);

                double d2Udn2 = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (dx * dx) +
                                (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dy * dy);
                double d2Vdn2 = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (dx * dx) +
                                (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dy * dy);

                U_new[i][j] = U[i][j] - dt * (U[i][j] * dUdn + V[i][j] * dVdn) + (dt / Re) * d2Udn2;
                V_new[i][j] = V[i][j] - dt * (U[i][j] * dVdn + V[i][j] * dVdn) + (dt / Re) * d2Vdn2;
        
            
            }
        }


        // Apply boundary conditions
            if (coords[1] == 0) { // x = 0
                for (int j = 1; j <= local_Ny; ++j) {
                    U[0][j] = 0.0;
                    U_new[0][j] = 0.0;

                    V[0][j] = 0.0;
                    V_new[0][j] = 0.0;
                }
            }

            if (coords[1] == dims[0] - 1) { // x = 1
                for (int j = 1; j <= local_Ny; ++j) {
                    double y = (global_y_start + j - 1) * dy;
                    if (y > 0.7 && y < 1.0) {
                        U[local_Nx + 1][j] = -1.0;
                        U_new[local_Nx + 1][j] = -1.0;

                        V[local_Nx + 1][j] = 0.0;
                        V_new[local_Nx + 1][j] = 0.0;
                    } else {
                        U[local_Nx + 1][j] = 0.0;
                        U_new[local_Nx + 1][j] = 0.0;

                        V[local_Nx + 1][j] = 0.0;
                        V_new[local_Nx + 1][j] = 0.0;
                    }
                }
            }

            if (coords[0] == 0) { // y = 0
                for (int i = 1; i <= local_Nx; ++i) {
                    double x = (global_x_start + i - 1) * dx;
                    if (x < 0.3) {
                        U[i][0] = 0.0;
                        U_new[i][0] = 0.0;

                        V[i][0] = 0.0;
                        V_new[i][0] = 0.0;
                    } else if (x >= 0.3 && x <= 0.6) {
                        U[i][0] = U[i][1]; // Neumann
                        U_new[i][0] = U_new[i][1];

                        V[i][0] = V[i][1];
                        V_new[i][0] = V_new[i][1];
                    } else {
                        U[i][0] = 0.0;
                        U_new[i][0] = 0.0;
                        V[i][0] = 0.0;
                        V_new[i][0] = 0.0;
                    }
                }
            }

            if (coords[0] == dims[1] - 1) { // y = 1
                for (int i = 1; i <= local_Nx; ++i) {
                    double x = (global_x_start + i - 1) * dx;
                    if (x < 0.3) {
                        U[i][local_Ny + 1] = U[i][local_Ny]; // Neumann
                        U_new[i][local_Ny + 1] = U_new[i][local_Ny];

                        V[i][local_Ny + 1] = V[i][local_Ny];
                        V_new[i][local_Ny + 1] = V_new[i][local_Ny];
                    } else {
                        U[i][local_Ny + 1] = 0.0;
                        U_new[i][local_Ny + 1] = 0.0;

                        V[i][local_Ny + 1] = 0.0;
                        V_new[i][local_Ny + 1] = 0.0;
                    }
                }
            }

       

        U = U_new;
        V = V_new;
    }


    save_results(U, V, local_Nx, local_Ny, global_x_start, global_y_start, rank, size, dims, cart_comm);
    MPI_Finalize();
    return 0;
}
