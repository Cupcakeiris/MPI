#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

void save_to_file(const vector<float> &arr) {
    ofstream outfile("analytical_x.txt");
    for (const auto& val : arr) {
        outfile << val << "\n";
    }
    outfile.close();
}

int main() {
    int Nx = 1000;
    float L = 1.0;
    vector<float> analytical_x(Nx, 0);

    for (int i = 0; i < Nx; i++) {
        float dx = L / (Nx - 1); 
        float x = i * dx;
        analytical_x[i] = -exp(1) * cos(1) * x + exp(x) * cos(x);
    }

    save_to_file(analytical_x);
    cout << "Data saved to analytical_x.txt" << endl;

    return 0;
}

