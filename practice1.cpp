#include <iostream>
#include <vector>
using namespace std;

//Ready formulas
float S(float a, float p)  {return 1/((1-a) + a/p);}
float acceleration(float a, float b){return a/b;}
float timeSpent(float a, float b){return a*b;}
float efficiency(float a, float b, float p){return a/b/p;}
float product(vector<vector<float>>& A, vector<vector<float>>& B){
    int rowA = A.size();
    int colA = A[0].size();
    int rowB = B.size();
    int colB = B[0].size();

    if(colB != rowA){printf("Wrong dimensions"); return 0;}
    float res[rowA][colB];
    for(int i = 0; i < rowA; i++){
        for(int j = 0; j < colB; j++){
            res[i][j] = 0;          //All elements will be 0, to perform dot product of vector rowA and vector colB, they have k elements
            for(int k=0; k<rowB; k++){
                res[i][j] += A[i][k] * B[i][k];
            }
            printf("%f\t", res[i][j]);
        }
        printf("\n");
    }
    return 0;
}

int main(){
    // FLOPS
    float flops, tflops; // transferred flops
    flops = 50e9;  tflops = flops/1e12; printf("%f TFLOPS\n", tflops);
    flops = 50e12; tflops = flops/1e9;  printf("%f GFLOPS\n", tflops);
    flops = 17e15; tflops = flops/1e15; printf("%f PFLOPS\n", tflops);
    flops = 15e12; tflops = flops/1e15; printf("%f PFLOPS\n", tflops);
    flops = 25e12; tflops = flops/1e6;  printf("%6.0f MFLOPS\n", tflops);

    cout<<endl<<endl;

    //Amdahl's law: S = 1/(a + (1-a)/p)
    printf("%f\n%f\n%f", 
            S(0.75, 4), 
            S(0.9, 10), 
            S(0.8, 6));

    cout<<endl<<endl;

    //Acceleration
    printf("%f\n%f\n%f", 
            acceleration(10, 3), 
            timeSpent(20, 4), 
            efficiency(15, 3, 5));

    cout<<endl<<endl;

    //Matrices multiplication
    //Give input inside the code rather from terminal, it's much optimized way (c) T. Mustafin while giving his first lab
    vector<vector<float>> A = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9},
        {10, 11, 12}
    };

    vector<vector<float>> B = {
        {0.5, 0, 0, 0},
        {0, 0.5, 0, 0},
        {0, 0, 0.5, 0},
    };

    product(A, B);
}