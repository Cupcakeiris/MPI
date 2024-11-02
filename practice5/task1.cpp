#include <iostream>
#include <mpi.h>
#include <cmath>

using namespace std;
double INTEGRAND(double x);
double factorial(int N);
double p_piii(int N, int rank, int size);
double piii(int N);
double p_exppp(int N, int rank, int size);
double exppp(int N);
double p_integgg(double lower, double upper, int N, int rank, int size);
double integgg(double lower, double upper, int N);



int main(int argc, char * argv[]){
// ///-------------------------------------------------------- Sequential program
    // // int N = 1e3;
    // double N = 1e5;
    // double res, start, end;
    // int a = 0, b = 10;

    // start = MPI_Wtime();
    // res = piii(N);
    // end = MPI_Wtime();
    // cout<<"Result: "<<res<<"\tTime taken: "<<end-start<<"\n";

    // start = MPI_Wtime();
    // res = exppp(N);
    // end = MPI_Wtime();
    // cout<<"Result: "<<res<<"\tTime taken: "<<end-start<<"\n";

    // start = MPI_Wtime();
    // res = integgg(a, b, N);
    // end = MPI_Wtime();
    // cout<<"Result: "<<res<<"\tTime taken: "<<end-start<<"\n";

// --------------------------------------------------------- Parallel program
    // int N = 1e3;
    double N = 1e5;
    double res_exp, res_pi, res_integ, start, end;
    double a = 0, b = 10;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Wtime();
    MPI_Request request;

    double global_res_pi = 0, global_res_exp = 0, global_res_integ = 0;
    double local_res_pi, local_res_exp, local_res_integ;

    res_pi = p_piii(N, rank, size);
    res_exp = p_exppp(N, rank, size);
    res_integ = p_integgg(a, b, N, rank, size);


//----------------------------------------------- Blocking communication
    // if(rank==0){ 
    //     // pi
    //     global_res_pi = res_pi;
    //     start = MPI_Wtime();
    //     for(int i = 1; i<size; i++){
    //         MPI_Recv(&local_res_pi, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 31415, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         global_res_pi += local_res_pi;
    //     }
    //     end = MPI_Wtime();
    //     cout << "Result: " << global_res_pi * 4 << "\tTime taken: "<<end - start<<"\n";
    
    //     // exp
    //     global_res_exp = res_exp;
    //     start = MPI_Wtime();
    //     for(int i = 1; i<size; i++){
    //         MPI_Recv(&local_res_exp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 27182, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         global_res_exp += local_res_exp;
    //     }
    //     end = MPI_Wtime();
    //     cout << "Result: " << global_res_exp << "\tTime taken: "<<end - start<<"\n";
        
    //     // integ (sending with batches)
    //     global_res_integ = res_integ;
    //     start = MPI_Wtime();
    //     for(int i = 1; i<size; i++){
    //         MPI_Recv(&local_res_integ, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         global_res_integ += local_res_integ;
    //     }
    //     end = MPI_Wtime();
    //     cout << "Result: " << global_res_integ << "\tTime taken: "<<end - start<<"\n";
        
    
    // } else {
    //     MPI_Send(&res_pi, 1, MPI_DOUBLE, 0, 31415, MPI_COMM_WORLD);
    //     MPI_Send(&res_exp, 1, MPI_DOUBLE, 0, 27182, MPI_COMM_WORLD);
    //     MPI_Send(&res_integ, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    // }


// // //-------------------------------------------------- Non-blocking communication

    MPI_Request request_pi, request_exp, request_integ;
    if(rank==0){ 
        // pi
        global_res_pi = res_pi;
        start = MPI_Wtime();
        for(int i = 1; i<size; i++){
            MPI_Irecv(&local_res_pi, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 31415, MPI_COMM_WORLD, &request_pi);
            MPI_Wait(&request_pi, MPI_STATUS_IGNORE);
            global_res_pi += local_res_pi;
        }
        end = MPI_Wtime();
        cout << "Result: " << global_res_pi * 4 << "\tTime taken: "<<end - start<<"\n";
    
        // exp
        global_res_exp = res_exp;
        start = MPI_Wtime();
        for(int i = 1; i<size; i++){
            MPI_Irecv(&local_res_exp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 27182, MPI_COMM_WORLD, &request_exp);
            MPI_Wait(&request_exp, MPI_STATUS_IGNORE);
            global_res_exp += local_res_exp;
        }
        end = MPI_Wtime();
        cout << "Result: " << global_res_exp << "\tTime taken: "<<end - start<<"\n";
        
        // integ
        global_res_integ = res_integ;
        start = MPI_Wtime();
        for(int i = 1; i<size; i++){
            MPI_Irecv(&local_res_integ, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request_integ);
            MPI_Wait(&request_integ, MPI_STATUS_IGNORE);
            global_res_integ += local_res_integ;
        }
        end = MPI_Wtime();
        cout << "Result: " << global_res_integ << "\tTime taken: "<<end - start<<"\n";
        
    
    } else {
        MPI_Isend(&res_pi, 1, MPI_DOUBLE, 0, 31415, MPI_COMM_WORLD, &request_pi);
        MPI_Isend(&res_exp, 1, MPI_DOUBLE, 0, 27182, MPI_COMM_WORLD, &request_exp);
        MPI_Isend(&res_integ, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request_integ);
        
    }



    MPI_Finalize();
}

double INTEGRAND(double x){
    return pow(x,3) - 1;
}

double factorial(int N){
    if (N == 0 || N == 1) {
        return 1.0;
    }
    double res = 1;
    for (int i = 2; i <= N; i++) {
        res *= i;
    }
    return res;
}

double p_piii(int N, int rank, int size){
    double sum = 0;
    for(int i = rank; i < N; i += size){
        sum += pow(-1, i) / (2 * i + 1);  

    }

    return sum;
}

double piii(int N){
    double sum = 0;
    for(int i=0; i<N; i++){
        sum += pow(-1, i)/(2*i+1);
    }
    return sum*4;
}

double p_exppp(int N, int rank, int size) {
    double sum = 0;
    for (int i = rank; i < N; i += size) { 
        sum += 1.0 / factorial(i);
    }
    return sum;
}


double exppp(int N){
    double sum = 0;
    for(int i=0; i<N; i++){
        sum += 1/factorial(i);
    }
    return sum;
}

double integgg(double lower, double upper, int N){
    double h = (upper - lower)/N;
    double sum = 0;

    for(double i = lower; i<upper; i += h){     // Using middle point Riemann sum
        sum += INTEGRAND(i)*h;
    }
    
    return sum;
}

double p_integgg(double lower, double upper, int N, int rank, int size){
    double h = (upper - lower) / N;   
    double local_sum = 0.0;
    
    for (int i = rank; i < N; i += size) {
        double x = lower + (i + 0.5)*h; 
        local_sum += INTEGRAND(x)*h;   
        
    }

    return local_sum;
}