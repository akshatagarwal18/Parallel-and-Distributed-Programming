#include "mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

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
    } return;
}

void print_matrix(float* A, int row, int col){
    for(int i = 0 ; i < row ; i++){
        for(int j = 0 ; j < col ; j++){
            printf("%6.3f ",A[i*row+j]);
        }
    }
}

float* create_array(int n, int m){
    float* mat = (float*)malloc(n*m*sizeof(float)); return mat;
}

void print_error(int n , float* A, float* B, float* C){
    float* D = create_array(n,n);
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            for(int k = 0 ; k < 32 ; k++){
                D[i*n + j] = D[i*n+j] + A[i*32 + k]*B[k*n + j];
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
    float *A, *B, *C;
    int tag,row_partition, first, extra_rows, total_elems, send_rows;
    MPI_Init(&argc, &argv);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    if(size <= 1){
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(1);
    }

    if(rank == MASTER){
        A = create_array(n,THIRTY_TWO);
        B = create_array(THIRTY_TWO,n);
        C = create_array(n,n);
        srand((unsigned int)time(NULL));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < 32 ; j++){
                A[i*32 + j] = rand()/(float)RAND_MAX ;
            }
        }
        for(int i = 0 ; i < 32 ; i++){
            for(int j = 0 ; j < n ; j++){
                B[i*n + j] = rand()/(float)RAND_MAX ;
            }
        }

        row_partition = n/size;
        extra_rows = n%size;
        tag = MASTER_TO_SLAVE_TAG;
        total_elems = 0;
        double start = MPI_Wtime();
        for(int i = 0 ; i < size; i++){
            send_rows = row_partition;
            if(i < extra_rows){send_rows++;}
            if(i == 0){
                first = send_rows;
                total_elems+=send_rows;
                continue;
            }
            MPI_Send(&total_elems,1,MPI_INT,i,tag,MPI_COMM_WORLD);
            MPI_Send(&send_rows,1,MPI_INT,i,tag,MPI_COMM_WORLD);
            MPI_Send(A+THIRTY_TWO*total_elems, send_rows*THIRTY_TWO, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
            MPI_Send(B, THIRTY_TWO*n, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
            total_elems+=send_rows;
        }

        for(int i = 0 ; i < first ; i++){
            for(int j = 0 ; j < n ; j++){
                for(int k = 0 ; k < 32 ; k++){
                    C[i*n+j]+=(A[i*32+k]*B[k*n+j]);
                }
            }
        }

        tag = SLAVE_TO_MASTER_TAG;
        for(int i = 1 ; i < size ; i++){
            MPI_Recv(&total_elems,1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
            MPI_Recv(&send_rows,1,MPI_INT,i,tag,MPI_COMM_WORLD,&status);
            MPI_Recv(C+n*total_elems,send_rows*n,MPI_FLOAT,i,tag,MPI_COMM_WORLD,&status);
        }
        double end = MPI_Wtime();
        double tim = end-start;

        print_matrix(A,n,32);
        print_matrix(B,32,n);
        print_matrix(C,n,n);
        print_error(n,A,B,C);
        printf("Time taken for matrix multiplication is %0.4fs\n ", tim);

    }
    else{
        tag = MASTER_TO_SLAVE_TAG;
        MPI_Recv(&total_elems, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&send_rows, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD, &status);

        A = create_array(send_rows,THIRTY_TWO);
        B = create_array(THIRTY_TWO,n);
        C = create_array(send_rows,n);

        MPI_Recv(A, send_rows*THIRTY_TWO, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(B, THIRTY_TWO*n, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);

        for(int i = 0 ; i < send_rows ; i++){
            for(int j = 0 ; j < n ; j++){
                for(int k = 0 ; k < 32 ; k++){
                    C[i*n+j]+=(A[i*32+k]*B[k*n+j]);
                }
            }
        }
        tag = SLAVE_TO_MASTER_TAG;

        MPI_Send(&total_elems, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD);
        MPI_Send(&send_rows, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD);
        MPI_Send(C,send_rows*n, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}