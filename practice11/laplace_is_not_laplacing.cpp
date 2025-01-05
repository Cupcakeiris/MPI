#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

void save_to_file(const vector<vector<double>> &local_u, int rank, int nx, int ny, int coords[], double hx, double hy) {
    string filename = "laplace_" + to_string(rank) + ".txt";
    ofstream fout(filename);

    if (!fout.is_open()) {
        cerr << "Error opening fil: " << filename << endl;
        return;
    }

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double x = (coords[0] * nx + (i - 1)) * hx;
            double y = (coords[1] * ny + (j - 1)) * hy;
            fout << x << " " << y << " " << local_u[i][j] << "\n";
        }
    }

    fout.close();
}


const double Lx = 2.0, Ly = 1.0;
const int Nx = 100, Ny = 50;
const double tol = 1e-6;
const int max_iter = 10000;

double boundary_condition(double x, double y) {
    if (x == 2.0) return 2 * y * y + y;
    if (y == 1.0) return 0.5 * x * x + 0.5 * x;
    return 0.0;
}

void exchange_boundaries(vector<vector<double>> &local_u, int rank, int up, int down, int left, int right, MPI_Comm comm) {
    int nx = local_u.size() - 2;
    int ny = local_u[0].size() - 2;

    vector<double> send_up(ny), send_down(ny), recv_up(ny), recv_down(ny);
    vector<double> send_left(nx), send_right(nx), recv_left(nx), recv_right(nx);

    if (up != MPI_PROC_NULL) {
        for (int j = 0; j < ny; ++j)
            send_up[j] = local_u[1][j + 1];
    }
    if (down != MPI_PROC_NULL) {
        for (int j = 0; j < ny; ++j)
            send_down[j] = local_u[nx][j + 1];
    }
    if (left != MPI_PROC_NULL) {
        for (int i = 0; i < nx; ++i)
            send_left[i] = local_u[i + 1][1];
    }
    if (right != MPI_PROC_NULL) {
        for (int i = 0; i < nx; ++i)
            send_right[i] = local_u[i + 1][ny];
    }

    MPI_Sendrecv(send_up.data(), ny, MPI_DOUBLE, up, 0, recv_down.data(), ny, MPI_DOUBLE, down, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_down.data(), ny, MPI_DOUBLE, down, 1, recv_up.data(), ny, MPI_DOUBLE, up, 1, comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv(send_left.data(), nx, MPI_DOUBLE, left, 2, recv_right.data(), nx, MPI_DOUBLE, right, 2, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_right.data(), nx, MPI_DOUBLE, right, 3, recv_left.data(), nx, MPI_DOUBLE, left, 3, comm, MPI_STATUS_IGNORE);

    if (up != MPI_PROC_NULL) {
        for (int j = 0; j < ny; ++j)
            local_u[0][j + 1] = recv_up[j];
    }
    if (down != MPI_PROC_NULL) {
        for (int j = 0; j < ny; ++j)
            local_u[nx + 1][j + 1] = recv_down[j];
    }
    if (left != MPI_PROC_NULL) {
        for (int i = 0; i < nx; ++i)
            local_u[i + 1][0] = recv_left[i];
    }
    if (right != MPI_PROC_NULL) {
        for (int i = 0; i < nx; ++i)
            local_u[i + 1][ny + 1] = recv_right[i];
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm comm_cart;
    int dims[2] = {0, 0};
    int periods[2] = {0, 0};

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Cartesian topology
    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);

    int coords[2], up, down, left, right;
    MPI_Cart_coords(comm_cart, rank, 2, coords);
    MPI_Cart_shift(comm_cart, 0, 1, &up, &down);
    MPI_Cart_shift(comm_cart, 1, 1, &left, &right);

    // Size of local net
    int nx = Nx / dims[0];
    int ny = Ny / dims[1];
    double hx = Lx / Nx, hy = Ly / Ny;

    vector<vector<double>> local_u(nx + 2, vector<double>(ny + 2, 0.0));

    // Boundary conditions
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double x = (coords[0] * nx + (i - 1)) * hx;
            double y = (coords[1] * ny + (j - 1)) * hy;

            if (coords[0] == 0 && i == 1)
                local_u[i][j] = boundary_condition(0, y);
            if (coords[0] == dims[0] - 1 && i == nx)
                local_u[i + 1][j] = boundary_condition(Lx, y);
            if (coords[1] == 0 && j == 1)
                local_u[i][j] = boundary_condition(x, 0);
            if (coords[1] == dims[1] - 1 && j == ny)
                local_u[i][j + 1] = boundary_condition(x, Ly);
        }
    }

    // Gauss-Seidel
    double diff = tol + 1;
    int iter = 0;

    while (diff > tol && iter < max_iter) {
        iter++;
        exchange_boundaries(local_u, rank, up, down, left, right, comm_cart);

        diff = 0.0;
        for (int i = 1; i <= nx; ++i) {
            for (int j = 1; j <= ny; ++j) {
                double old_u = local_u[i][j];
                local_u[i][j] = 0.5 * ((local_u[i + 1][j] + local_u[i - 1][j]) / (hx * hx) +
                                       (local_u[i][j + 1] + local_u[i][j - 1]) / (hy * hy)) /
                                (1.0 / (hx * hx) + 1.0 / (hy * hy));

                diff = max(diff, fabs(local_u[i][j] - old_u));
            }
        }

        // To get max error
        double global_diff;
        MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
        diff = global_diff;
    }

    if (rank == 0) {
        cout << "Iterations to converge: " << iter << endl;
    }

    save_to_file(local_u, rank, nx, ny, coords, hx, hy);
    
    MPI_Barrier(comm_cart);

    MPI_Finalize();
    return 0;
}