#include "mpi.h"
#include<bits/stdc++.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
using namespace std;

#define MASTER 0
#define MASTER_TO_SLAVE_TAG 0
#define SLAVE_TO_MASTER_TAG 1
#define THIRTY_TWO 32

int isEqual(float* A , float* B, int row, int col){
    for(int i = 0 ; i < row ; i++){
        for(int j = 0 ; j < col ; j++){
            if(A[i*row + j]!=B[i*row + j]){return 0;}
        }
    } return 1;
}

void matrix_multiply(float* A, float* B, float* C, int m1, int n1, int p1){
    for(int i = 0 ; i < m1 ; i++){
        for(int j = 0 ; j < p1 ; j++){
            C[i*p1 + j] = 0.0;
            for(int k = 0 ; k < n1 ; k++){
                C[i*p1 + j]+=(A[i*n1+k]*B[k*p1+j]);
            }
        }
    } 

    return;
}

void print_matrix(float* A, int row, int col){
    for(int i = 0 ; i < row ; i++){
        for(int j = 0 ; j < col ; j++){
            //printf("%6.3f ",A[i*row+j]);
            cout<<A[i*row+j]<<" ";
        }cout<<endl;
    } cout<<endl;
}

float* create_array(int n, int m){
    float* mat = (float*)malloc(n*m*sizeof(float)); return mat;
}

void print_error(int n , float* A, float* B, float* C){
    float* D = create_array(n,n);
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
             D[i*n + j] = 0.0 ;
            for(int k = 0 ; k < THIRTY_TWO ; k++){
                D[i*n + j] = D[i*n+j] + A[i*THIRTY_TWO + k]*B[k*n + j];
            }
        }
    } float err = 0.0;
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            err+=fabs(D[i*n+j]-C[i*n+j]);
        }
    }
    printf("Error is %.14lf\n",err);
}

int main(int argc, char* argv[]){
    int n = 50;
    if(argc > 1){n=atoi(argv[1]);}
   // printf("n=%d\n",n);
    //float *A, *B, *C;
    int tag,row_partition , first, extra_rows, total_elems, send_rows;
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    row_partition = n/size;
    extra_rows = n%size;
    float *A, *B, *C;
    float *locA = create_array(row_partition,THIRTY_TWO);
    float *locC = create_array(row_partition,n);
    
    if(rank == MASTER){
        
        A = create_array(n,THIRTY_TWO);
        B = create_array(THIRTY_TWO,n);
        C = create_array(n,n);
        srand((unsigned int)time(NULL));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < THIRTY_TWO ; j++){
                A[i*THIRTY_TWO + j] = rand()/(float)RAND_MAX ;
            }
        }
        for(int i = 0 ; i < THIRTY_TWO ; i++){
            for(int j = 0 ; j < n ; j++){
                B[i*n + j] = rand()/(float)RAND_MAX ;
            }
        }
        double start = MPI_Wtime();
        MPI_Scatter(A,THIRTY_TWO*row_partition,MPI_FLOAT,locA,THIRTY_TWO*row_partition,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Bcast(B,THIRTY_TWO*n,MPI_FLOAT,0,MPI_COMM_WORLD);
        matrix_multiply(locA,B,locC,row_partition,THIRTY_TWO,n);
        MPI_Gather(locC,n*row_partition,MPI_FLOAT,C,n*row_partition,MPI_FLOAT,0,MPI_COMM_WORLD);
        double end = MPI_Wtime();
        double time = end-start;
        printf ("Time = %0.4fs \n",time);
        // print_matrix(A,n,THIRTY_TWO);
        // print_matrix(B,THIRTY_TWO,n);
        // print_matrix(C,n,n);
        print_error(n,A,B,C);
    }else{
        B = create_array(THIRTY_TWO,n);
        MPI_Scatter(A,THIRTY_TWO*row_partition,MPI_FLOAT,locA,THIRTY_TWO*row_partition,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Bcast(B,THIRTY_TWO*n,MPI_FLOAT,0,MPI_COMM_WORLD);
        matrix_multiply(locA,B,locC,row_partition,THIRTY_TWO,n);
        MPI_Gather(locC,n*row_partition,MPI_FLOAT,C,n*row_partition,MPI_FLOAT,0,MPI_COMM_WORLD);    
    }
    MPI_Barrier(MPI_COMM_WORLD);    
    MPI_Finalize();
}