#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
    const int nx = 100; 
    const int nt = 100; 
    const double L = 1.0; 
    const double T = 1.0; 
    const double dx = L / (nx - 1); 
    const double dt = T / nt; 
    const double c = 1.0; 

    vector<double> x(nx);
    vector<vector<double>> u(nt, vector<double>(nx));

    for (int i = 0; i < nx; ++i) {
        x[i] = i * dx;
    }

    // Initial condition: u_0(x)
    for (int i = 0; i < nx; ++i) {
        if (x[i] > 0) {
            u[0][i] = 1.0;
        } else {
            u[0][i] = 0.0;
        }
    }

    // Lax-Wendroff method
    for (int n = 0; n < nt - 1; ++n) {
        for (int i = 1; i < nx - 1; ++i) {
            u[n + 1][i] = u[n][i] - (c * dt / (2 * dx)) * (u[n][i + 1] - u[n][i - 1]) +
                          (c * c * dt * dt / (2 * dx * dx)) * (u[n][i + 1] - 2 * u[n][i] + u[n][i - 1]);
        }
        u[n + 1][0] = 0.0;
    }


    ofstream outfile("lax_wendroff_solution.txt");
    if (!outfile.is_open()) {
        cerr << "Error opening file" << endl;
        return 1;
    }


    for (int n = 0; n < nt; ++n) {
        double t = n * dt;
        for (int i = 0; i < nx; ++i) {
            outfile << t << " " << x[i] << " " << u[n][i] << endl;
        }
        outfile << endl;
    }

    outfile.close();
    cout << "Saved to lax_wendroff_solution.txt" << endl;

    return 0;
}