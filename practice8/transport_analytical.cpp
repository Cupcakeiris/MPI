#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {
    const int Nx = 100; 
    const int Nt = 100; 
    const double L = 1.0; 
    const double T = 1.0;
    const double dx = L / Nx;
    const double dt = T / Nt;

    vector<vector<double>> u(Nt, vector<double>(Nx));

    for (int n = 0; n < Nt; ++n) {
        double t = n * dt;
        for (int i = 0; i < Nx; ++i) {
            double current_x = i * dx;
            if (current_x > t) {
                u[n][i] = 1.0;
            } else {
                u[n][i] = 0.0;
            }
        }
    }

    ofstream outfile("transport_solution.txt");
    if (!outfile.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    for (int n = 0; n < Nt; ++n) {
        double t = n * dt;
        for (int i = 0; i < Nx; ++i) {
            outfile << t << " " << i * dx << " " << u[n][i] << endl;
        }
        outfile << endl; 
    }

    outfile.close();
    cout << "Saved to transport_solution.txt" << endl;

    return 0;
}