#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
using namespace std;

double pi = 3.14159265358979323846;

// Reading parameters from file
void read_file_input(const string &filename, double &L, double &T, double &c, double &Nx, double &Nt, double &t0) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line;

    while (getline(file, line)) {

        // cout << "Reading line: " << line << endl;

        istringstream iss(line);
        string parameter;
        char eq;
        double value;

        if (iss >> parameter >> eq >> value) {
            if (parameter == "L")       L = value;
            else if (parameter == "T")  T = value;
            else if (parameter == "c")  c = value;
            else if (parameter == "Nx") Nx = value;
            else if (parameter == "Nt") Nt = value;
            else if (parameter == "t0") t0 = value;
        }

    }
}

void save_to_file(const vector<double> &u, const string &filename) {
    ofstream outfile(filename);
    for (double val : u) {
        outfile << val << "\n";
    }
    outfile.close();
}

double initial_condition(double x) {
    return cos(pi * x);
}

double boundary_condition(double t) {
    return exp(-t);
}

void solve_backward_difference(double t0, int rank, int size, double dx, double dt, double alpha, int Nx) {
    int k = Nx / size;
    int i1 = rank * k;
    vector<double> u(Nx), u_new(Nx);
    vector<double> local_u(k + 2, 0), local_u_new(k + 2, 0);

    for (int i = 0; i < Nx; i++) {
        u[i] = initial_condition(i * dx);
    }

    MPI_Scatter(u.data(), k, MPI_DOUBLE, &local_u[1], k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double t = 0.0;
    while (t < t0) {
        if (rank == 0) {
            local_u[0] = boundary_condition(t);
        }

        if (rank < size - 1) {
            MPI_Send(&local_u[k], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&local_u[k + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank > 0) {
            MPI_Send(&local_u[1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&local_u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i <= k; i++) {
            local_u_new[i] = local_u[i] - alpha * (local_u[i] - local_u[i - 1]);
        }

        local_u = local_u_new;
        t += dt;
    }

    MPI_Gather(&local_u[1], k, MPI_DOUBLE, u.data(), k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        save_to_file(u, "solution_backward.txt");
        cout << "Saved to solution_backward.txt\n";
    }
}

void solve_analytical(double t0, int rank, int size, double dx, int Nx) {
    int k = Nx / size;
    int i1 = rank * k;
    vector<double> u(Nx);

    for (int i = 0; i < k; i++) {
        double x = (i1 + i) * dx;
        u[i] = cos(pi * (x - t0)) + exp(-t0) * (1 - x);
    }
   
    if (rank == 0) {
        save_to_file(u, "analytical.txt");
        cout << "Saved to analytical.txt\n";
    }
}

int main(int argc, char** argv) {
    double L, T, c, Nx, Nt, t0;
    read_file_input("params.txt", L, T, c, Nx, Nt, t0);

    double dx = L / (Nx - 1);
    double dt = T / Nt;
    double alpha = c * dt / dx;
    bool stability_check = (alpha <= 1.0);

    double start;
    double end;
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (!stability_check && rank == 0) {
        cout << "Stability condition alpha <= 1.0 is not met\n";
    }

    solve_analytical(t0, rank, size, dx, Nx);

    if(rank == 0){
        start = MPI_Wtime();
    }
    solve_backward_difference(t0, rank, size, dx, dt, alpha, Nx);
    
    if(rank == 0){
        end = MPI_Wtime();
        cout<<"Time taken: "<<end - start<<endl;
    }

    MPI_Finalize();
    return 0;
}
