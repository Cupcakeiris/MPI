#include <mpi.h>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

double analytical();
double phi_1(double x);
double phi_2(double x);
double f(double x, double y);


double double_Riemann(int rank, int size, double a, double b, int M, int N) {
    double dx = (b - a) / (M - 1);
    double dy;
    double total = 0;

    for(int i=rank; i < M - 1; i+=size){
        double xi = a + i*dx;
        double bb = phi_2(xi);
        double aa = phi_1(xi);
        dy = (bb - aa)/(N-1);
        double area = dx * dy;
        
        for(int j = 0; j < N - 1; j ++){
            double yj = aa + j*dy;
            total += f(xi, yj) * area;
        }
    }

    return total;
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a = 1;
    double b = 2;
    int M = 1000;
    int N = 1000;
    double numerical;
    double start, end;

    MPI_Bcast(&N, 1, MPI_INT, 0, comm);
    MPI_Bcast(&M, 1, MPI_INT, 0, comm);

    double partial_sum = double_Riemann(rank, size, a, b, N, M);  // Each process computes its part

    vector<double> all_partial_sums(size);                     // To store all partial sums
    MPI_Gather(&partial_sum, 1, MPI_DOUBLE, all_partial_sums.data(), 1, MPI_DOUBLE, 0, comm);  // Gather partial sums

    if (rank == 0) {
        start = MPI_Wtime();
        numerical = 0.0;
        for (int i = 0; i < size; i++) {
            numerical += all_partial_sums[i];
        }
        end = MPI_Wtime();
        double exact_value = analytical();
        double error = exact_value - numerical;
        cout << "Numerical: " << numerical << "\n";
        cout << "Exact: " << exact_value << "\n";
        cout << "Error: " << error << "\n";
        cout << "Time taken: " << end - start;
    }

    MPI_Finalize();
    return 0;
}

double analytical() {
    return 4 * log(4) - log(2) - 4;
}

double phi_1(double x) {
    return 1 - x;
}

double phi_2(double x) {
    return x;
}

double f(double x, double y) {
    return log(x + y) / x;
}