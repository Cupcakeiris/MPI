#include <mpi.h>
#include <iostream>
#include <vector>
using namespace std;

void print_matrix(const vector<vector<int>>& A);
int product(const vector<vector<int>>& A, const vector<int>& B) ;
void print_vector(const vector<int>& vec);

int main(int argc, char * argv[]){
    int rank, size, error_code;
    int number_of_working_ranks = 3;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Wtime(); // Time starts when MPI is initialized
    MPI_Status status;



    ///////////////////////////////////////////////////// Task 1
    // int ping;
    // if(rank==0){
    //     ping = 1337;
    //     MPI_Ssend(&ping, 1, MPI_INT, 1, 11, MPI_COMM_WORLD);
    //     cout<<"Process 0 Ssent "<<ping<<" to Process 1\n";
    //     double start_ssend = MPI_Wtime(); // Suppose time is now 0.01 sec from initialization

    //     MPI_Bsend(&ping, 1, MPI_INT, 2, 12, MPI_COMM_WORLD);
    //     cout<<"Process 0 Bsent "<<ping<<" to Process 2\n";
    //     double start_bsend = MPI_Wtime();

    //     MPI_Rsend(&ping, 1, MPI_INT, 3, 13, MPI_COMM_WORLD);
    //     cout<<"Process 0 Rsent "<<ping<<" to Process 3\n";
    //     double start_rsend = MPI_Wtime();
        
    // }
    // if(rank==1 || rank==2 || rank==3){
    //     error_code = MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    //     if(error_code != MPI_SUCCESS){cout<<"Process "<< rank <<" wasn't probed due to error\t"<<error_code<<"\n";}
    //     else{
    //         cout<<"Process "<<rank<<" was probed\t"<<"\n";
    //         MPI_Recv(&ping, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         cout<<"Process "<<rank<<" received "<<ping<<" from Process 0\t"<<"| Time elapsed: "<<MPI_Wtime()<<"\n";
    //     // and that now it's 0.03 sec from initialization
    //     // result of receiving after sending is 0.03s - 0.01s = 0.02s
    //     }
        
    // }




    ///////////////////////////////////////////////////// Task 2
    // int sum;
    // if(rank==0){
    //     vector<int> pong = {1, 3, 3, 7};
    //     vector<int> peng = {1, 3, 3, 3, 3, 7};
    //     vector<int> pang = {1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7};

    //     MPI_Send(pong.data(), pong.size(), MPI_INT, 1, 0, MPI_COMM_WORLD);
    //     MPI_Send(peng.data(), peng.size(), MPI_INT, 2, 0, MPI_COMM_WORLD);
    //     MPI_Send(pang.data(), pang.size(), MPI_INT, 3, 0, MPI_COMM_WORLD);

    //     MPI_Recv(&sum, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     cout<<"Process 0 received "<<sum<<" from process "<<1<<"\n";
    //     MPI_Recv(&sum, 1, MPI_INT, 2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     cout<<"Process 0 received "<<sum<<" from process "<<2<<"\n";
    //     MPI_Recv(&sum, 1, MPI_INT, 3, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     cout<<"Process 0 received "<<sum<<" from process "<<3<<"\n";
    // }    

    // if(rank==1 || rank ==2 || rank==3){
    //     MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    //     int count;
    //     MPI_Get_count(&status, MPI_INT, &count);
    //     int recv[count]; // We calculate how many space we need by probing
    //     MPI_Recv(recv, count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //     sum = 0;
    //     for(int i=0; i<count; i++){sum+=recv[i];}
    //     MPI_Send(&sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    //     cout<<"Process "<< rank << " sent sum: "<< sum<<"\n";
    // }
    


    ///////////////////////////////////////////////////// Task 3
    int n = 7;
    int m = 6;
    vector<vector<int>> A(n, vector<int>(m)); // Ax = b
    vector<int> x(m);
    vector<int> b(n);

    if(rank==0){
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = (rand() % 10) + 1;
            }
        }
        for (int i = 0; i < m; i++) {
            x[i] = (rand() % 10) + 1;
        }
                
        cout<<"-----------Generated values:---------------\n\n";
        print_matrix(A);
        cout<<"\n";
        print_vector(x);

        // sends each row
        for (int worker_rank = 1; worker_rank < size; worker_rank++) {
            MPI_Send(A[worker_rank - 1].data(), A[0].size(), MPI_INT, worker_rank, 0, MPI_COMM_WORLD);  
            MPI_Send(x.data(), A[0].size(), MPI_INT, worker_rank, 1, MPI_COMM_WORLD); 
        }

        double start_time = MPI_Wtime();
        for (int worker_rank = 1; worker_rank < size; worker_rank++) {
            MPI_Recv(&b[worker_rank - 1], 1, MPI_FLOAT, worker_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        cout<<"------------------Product result--------------\n\n";
        double end_time = MPI_Wtime();
        print_vector(b);
        cout<<"Time spent: "<<end_time - start_time<<"\n";      
        cout<<"--------------------Analytical result-----------\n\n";
        double start_analytical = MPI_Wtime();
        product(A, x);
        double end_analytical = MPI_Wtime();
        cout<<"Time spent: "<<end_analytical - start_analytical<<"\n";
    }

    else{
        vector<int> A_row(A[0].size()); // in 1 row we have # of columns elements
        vector<int> x_copy(x.size());

        MPI_Recv(A_row.data(), A[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(x_copy.data(), A[0].size(), MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int result = 0;
        for (int i = 0; i < A[0].size(); i++) {
            result += A_row[i] * x_copy[i];
        }
        MPI_Send(&result, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }



    ///////////////////////////////////////////////////// Task 4
    // int pung, pyng;
    // MPI_Request send_request, recv_request;
    // MPI_Status stat;

    // if(rank==0){
    //     pung = 1337;

    // //// --------------------- Deadlock ------------------------
    //     // MPI_Send(&pung, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); 
    //     // cout<<"Process 0 sent "<<pung<<" to Process 1 \n";
    //     // MPI_Recv(&pyng, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     // cout<<"Process 0 received "<<pyng<<" from Process 1\n";
    
    // //// --------------------- Non deadlock ------------------------
    //     MPI_Isend(&pung, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &send_request); 
    //     cout<<"Process 0 sent "<<pung<<" to Process 1 \n";
    //     double start_time_0 = MPI_Wtime();

    //     MPI_Irecv(&pyng, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);
    //     cout<<"Process 0 is going to receive from Process 1\n";
    //     MPI_Wait(&recv_request, &stat);

    //     double end_time_0 = MPI_Wtime();
    //     cout<<"Process 0 received "<<pyng<<" from Process 1\t\t" << end_time_0 - start_time_0<<"\n";
    // }
    
    // if(rank==1){
    //     pyng = 3173;

    // //// --------------------- Deadlock ------------------------
    //     // MPI_Send(&pyng, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    //     // cout<<"Process 1 sent "<<pyng<<" to Process 0 \n";
    //     // MPI_Recv(&pung, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     // cout<<"Process 1 received "<<pung<<" from Process 0\n";
    
    // //// --------------------- Non deadlock ------------------------
    //     MPI_Isend(&pyng, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &send_request);
    //     cout<<"Process 1 sent "<<pyng<<" to Process 0 \n";
    //     double start_time_1 = MPI_Wtime();

    //     MPI_Irecv(&pung, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);
    //     cout<<"Process 1 is going to receive from Process 0\n";
    //     MPI_Wait(&recv_request, &stat);

    //     double end_time_1 = MPI_Wtime();
    //     cout<<"Process 1 received "<<pung<<" from Process 0\t\t" << end_time_1 - start_time_1<<"\n";
    
    // }

    // MPI_Finalize();
    // return 0;
}

void print_matrix(const vector<vector<int>>& A) {
    for (const auto& row : A) {
            for (int value : row) {
                cout << value << "\t";
            }
            cout << "\n";
        }
}

void print_vector(const vector<int>& vec){
    for (int value : vec) {
        cout << value << "\t";
    }
    cout << endl;
}

int product(const vector<vector<int>>& A, const vector<int>& B) {
    int rows_A = A.size(); 
    int cols_A = A[0].size();   
    int size_B = B.size();     

    vector<int> result(rows_A, 0);

    for (int i = 0; i < rows_A; i++) {
        result[i] = 0;
        for (int k = 0; k < cols_A; k++) {
            result[i] += A[i][k] * B[k];
        }
        cout << result[i] << "\t";
    }
    cout<<"\n";
    return 0;
}