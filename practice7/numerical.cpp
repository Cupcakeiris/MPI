#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <mpi.h>
using namespace std;

void save_to_file(const vector<float> &arr) {
    ofstream outfile("newP.txt");
    for (const auto& val : arr) {
        outfile << val << "\n";
    }
    outfile.close();
}

float f(float x) {
    return 2 * exp(x) * cos(x);
}

int main() {

    double start = MPI_Wtime();

    int Nx = 1000;
    float L = 1.0;
    vector<float> newP(Nx, 0);
    vector<float> oldP(Nx, 0);
    float dx = L / (Nx - 1);
    oldP[0] = 1;  
    oldP[Nx - 1] = 0;

    float tol = 1e-6; 
    float max_err = 100.0;

    while (max_err > tol) { 
        max_err = 0.0;       
        newP[0] = 1;
        newP[Nx - 1] = 0;
        
        for (int i = 1; i < Nx - 1; i++) {
            float x = i * dx;
            newP[i] = 0.5 * (oldP[i - 1] + oldP[i + 1] + dx * dx * f(x));
            
            max_err = max(max_err, abs(newP[i] - oldP[i]));
        }

        oldP = newP;
    }

    double end = MPI_Wtime();

    save_to_file(newP);
    cout << "Data saved to newP.txt" << endl;
    cout<< "Time taken: "<<end - start << endl;

    return 0;
}
