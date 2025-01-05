#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {
    const int Nx = 20;
    const int Ny = 20;
    const double cx = 0.5;      // Speed in x direction
    const double cy = 1.0;      // Speed in y direction
    const int Nt = (Nx * Ny)/(cx + cy);
    const double Lx = 1.0;
    const double Ly = 1.0;
    const double T = 1.0;
    const double dx = Lx / Nx;
    const double dy = Ly / Ny;
    const double dt = T / Nt;

    vector<vector<vector<double>>> u(Nt, vector<vector<double>>(Nx, vector<double>(Ny, 0.0)));

    // Initial condition
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = i * dx;
            double y = j * dy;
            u[0][i][j] = (x > 0 && y > 0) ? 1.0 : 0.0;
        }
    }

    // Finite difference scheme
    for (int n = 0; n < Nt - 1; ++n) {
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double ux = (u[n][i][j] - u[n][i - 1][j]) / dx;
                double uy = (u[n][i][j] - u[n][i][j - 1]) / dy;
                u[n + 1][i][j] = u[n][i][j] - dt * (cx * ux + cy * uy);
            }
        }
    }


    ofstream outfile("transport_solution_2d.txt");
    if (!outfile.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    for (int n = 0; n < Nt; ++n) {
        double t = n * dt;
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                double x = i * dx;
                double y = j * dy;
                outfile << t << " " << x << " " << y << " " << u[n][i][j] << endl;
            }
            outfile << endl;
        }
        outfile << endl;
    }

    outfile.close();
    cout << "Saved to transport_solution_2d.txt" << endl;

    return 0;
}
