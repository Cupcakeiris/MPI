#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

int main() {
    const int Nx = 1000;  
    const int Nt = 2000;
    const double L = 10.0; 
    const double T = 6.0; 
    const double dx = L / (Nx - 1); 
    const double c = 2.0; 
    const double dt = T / Nt;
    const double stable = dx / c;
    const double alpha = c * dt / dx;
    
    if (dt > stable) {
        cout << "Scheme is not stable\n";
        return 1;
    }

    vector<vector<double>> u(Nt, vector<double>(Nx));

    for (int n = 0; n < Nt; ++n) {
        double t = n * dt;
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            if (x > c*t) {
                u[n][i] = cos(x - c*t) * cos(x - c*t);
            } else {
                u[n][i] = 1.0;
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