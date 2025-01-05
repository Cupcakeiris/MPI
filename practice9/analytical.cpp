#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {
    const int Nx = 20;
    const int Ny = 20;
    const double cx = 1.0;
    const double cy = 1.0;
    const double Lx = 1.0;
    const double Ly = 1.0;
    const double T = 1.0;
    const int Nt = ((Nx * Ny) / (cx + cy));

    const double dx = Lx / (Nx - 1);
    const double dy = Ly / (Ny - 1);
    const double dt = T / (Nt - 1);

    ofstream outfile("transport_analytical_2D.txt");
    if (!outfile.is_open()) {
        cerr << "Error opening file" << endl;
        return 1;
    }


    for (int n = 0; n < Nt; ++n) {
        double t = n * dt;
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            for (int j = 0; j < Ny; ++j) {
                double y = j * dy;
                double u = (x > cx * t && y > cy * t) ? 1.0 : 0.0;
                outfile << t << " " << x << " " << y << " " << u << endl;
                // cout<<"u: "<<u<<"\tx: "<<x<<"\tcx * t: "<< cx * t<<"\ty: "<<y<<"\tcy*t: "<<cy*t<<endl;
            }
            outfile << endl;
        }
        outfile << endl;
    }

    outfile.close();
    cout << "Saved to transport_analytical_2D.txt" << endl;
    return 0;
}
