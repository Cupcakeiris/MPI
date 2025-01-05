#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

const double Lx = 1.0, Ly = 1.0;
const int Nx = 60, Ny = 60;
const double dx = Lx / Nx, dy = Ly / Ny;
const double t = 1.0;
const int Nt = 3000;
const double dt = t/Nt;
const double Re = 40.0;

double initial_u(double x, double y) { return 0.0; }
double initial_v(double x, double y) { return 0.0; }

double norm(double x, double y) {
    return sqrt(x * x + y * y);
}

void apply_boundary_conditions(vector<vector<double>> &U, vector<vector<double>> &V) {
    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {
            double x = i * dx, y = j * dy;

            // x = 1
            if (x == 1.0) {
                if (y > 0.7 && y < 1.0) {
                    U[i][j] = -1.0; // Dirichlet
                    V[i][j] = 0.0;
                } else {
                    U[i][j] = 0.0;
                    V[i][j] = 0.0; // Dirichlet
                }
            }

            // y = 1
            if (y == 1.0) {
                if (x >= 0 && x <= 0.3) {
                    if (j == Ny) {
                        U[i][j] = U[i][j - 1];
                        V[i][j] = V[i][j - 1];
                    }
                } else if (x > 0.3 && x <= 1.0) {
                    U[i][j] = 0.0;
                    V[i][j] = 0.0;
                }
            }

            // y = 0
            if (y == 0.0) {
                if (x >= 0 && x <= 0.3) {
                    U[i][j] = 0.0;
                    V[i][j] = 0.0;
                } else if (x > 0.3 && x <= 0.6) {
                    if (j == 0) { 
                        U[i][j] = U[i][j + 1];
                        V[i][j] = V[i][j + 1];
                    }
                } else if (x > 0.6 && x <= 1.0) {
                    U[i][j] = 0.0;
                    V[i][j] = 0.0;
                }
            }

            // x = 0
            if (x == 0.0) {
                U[i][j] = 0.0;
                V[i][j] = 0.0;
            }
        }
    }
}


int main() {
    vector<vector<double>> U(Nx + 1, vector<double>(Ny + 1, 0.0));
    vector<vector<double>> V(Nx + 1, vector<double>(Ny + 1, 0.0));
    vector<vector<double>> U_new = U, V_new = V;

    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {
            double x = i * dx, y = j * dy;
            U[i][j] = initial_u(x, y);
            V[i][j] = initial_v(x, y);
        }
    }

    for (int n = 0; n < Nt; ++n) {
        for (int i = 1; i < Nx; ++i) {
            for (int j = 1; j < Ny; ++j) {
                double x = i * dx, y = j * dy;
                double norm_n = norm(x, y);

                double dUdn = (U[i + 1][j] - U[i - 1][j]) * (0.5 * x / norm_n) +
                              (U[i][j + 1] - U[i][j - 1]) * (0.5 * y / norm_n);
                double dVdn = (V[i + 1][j] - V[i - 1][j]) * (0.5 * x / norm_n) +
                              (V[i][j + 1] - V[i][j - 1]) * (0.5 * y / norm_n);

                double d2Udn2 = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (dx * dx) +
                                (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dy * dy);
                double d2Vdn2 = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (dx * dx) +
                                (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dy * dy);

                // Explicit time-stepping
                U_new[i][j] = U[i][j] - dt * (U[i][j] * dUdn + V[i][j] * dVdn) + (dt / Re) * d2Udn2;
                V_new[i][j] = V[i][j] - dt * (U[i][j] * dVdn + V[i][j] * dVdn) + (dt / Re) * d2Vdn2;
            }
        }

        apply_boundary_conditions(U_new, V_new);

        U = U_new;
        V = V_new;
    }

    ofstream file("burgers_solution_norm.txt");
    for (int i = 0; i <= Nx; ++i) {
        for (int j = 0; j <= Ny; ++j) {
            double x = i * dx, y = j * dy;
            file << x << " " << y << " " << U[i][j] << " " << V[i][j] << "\n";
        }
        file << "\n";
    }
    file.close();

    cout << "Results saved to 'burgers_solution_norm.txt'" << endl;
    return 0;
}
