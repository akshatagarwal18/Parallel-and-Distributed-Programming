#include<bits/stdc++.h>
#include<omp.h>
#include<fstream>
#include <cstdlib>
#include <time.h>
#define endl "\n"
using namespace std;
using namespace chrono;


int n; // size of the matrix
int num_thread; //number of parallel threads
string testcase="";
bool randomize = false;
ofstream output_P, output_L, output_U;
double delta = 0.0000000001;

string print_matrix(double** matrix){
    ostringstream s;
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            s<<matrix[i][j]<<" ";
        }s<<endl;
    } return s.str();
}

double** create_matrix(){
    //double** matrix = new double*[n];
    double **matrix = (double**)malloc(n * sizeof(double*));
    for(int i = 0 ; i < n ; i++){
       // matrix[i] = new double[n];
       matrix[i] = (double*)malloc(n * sizeof(double));
        for(int j = 0 ; j < n ; j++){
            matrix[i][j]=0.0;
        }
    } return matrix;
}

void freeMemory(double** matrix){
    for(int i = 0 ; i < n ; i++){
        delete [] matrix[i];
    } delete [] matrix;
}

double cal_residue(double** P , double** A , double** L, double** U){
    double res = 0.0;
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            double temp = 0.0;
            for(int k = 0 ; k < n ; k++){
                double temp1 = P[i][k]*A[k][j];
                double temp2 = L[i][k]*U[k][j];
                temp = temp + (temp1-temp2);
            } res = res + temp*temp; 
        }
    } return res;
}

void initialize_A(double** A, double** AA){
    ifstream input(testcase);
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            if(randomize){
                A[i][j] = drand48()*100;
            }else{
                input>>A[i][j];
            }
            AA[i][j] = A[i][j];
        }
    }
}

void initialize_L(double** L){
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            if(i>j){L[i][j]=drand48()*100;}
            else if(i==j){L[i][j]=1.0;}
        }
    }
}

void initialize_U(double** U){
    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            if(i<=j){U[i][j]=drand48()*100;}            
        }
    }
}

void lu_decomposition(){
    output_P.open("P_MATRIX");
    output_L.open("L_MATRIX");
    output_U.open("U_MATRIX");
    // int* pi = new int(n);
    int* pi = (int*)malloc(n * sizeof(int));
    double **P = create_matrix();
    double **A = create_matrix(); 
    double **AA = create_matrix(); initialize_A(A,AA);
    double **L = create_matrix();  initialize_L(L);
    double **U = create_matrix();  initialize_U(U);

    auto t1 = high_resolution_clock::now();
    for(int i = 0 ; i < n ; i++){pi[i]=i;}

    for(int k = 0 ; k < n ; k++){
        double mx = 0.0;
        int k1 = -1;
        for(int i = k ; i < n ; i++){
            if(mx < abs(A[i][k])){
                mx = abs(A[i][k]);
                k1=i;
            }
        }
        if(abs(mx) < delta){
            perror("Received singular matrix A");
        }
        swap(pi[k],pi[k1]);
        swap(A[k],A[k1]);
        for(int i = 0 ; i < k ; i++){
            swap(L[k][i],L[k1][i]);
        }
        U[k][k]=A[k][k];
        for(int i = k+1 ; i < n ; i++){
            L[i][k] = A[i][k]/(U[k][k]);
            U[k][i] = A[k][i];
        }

        int i;
        #pragma omp parallel for num_threads(num_thread) private(i) shared(A,L,U)
        
            for(i = k+1 ; i<n ; i++)
            {
                for(int j = k+1; j<n ; j++)
                {
                    A[i][j] = A[i][j]-L[i][k]*U[k][j];
                }
            }
        
    }
    
        for(int j = 0 ; j < n ; j++){
            P[j][pi[j]]=1.0;
        }
        auto t2 = high_resolution_clock::now();
        cout<<"Time taken for LU Decomposition : "<<
        duration_cast<microseconds>(t2-t1).count()<<"\xC2\xB5s\n";
        cout << "Residue: " << cal_residue(P, AA, L, U)<<"\n";
        output_P<<print_matrix(P);
        output_L<<print_matrix(L);
        output_U<<print_matrix(U);
       

        freeMemory(P);        
        freeMemory(A);        
        freeMemory(AA);        
        freeMemory(L);        
        freeMemory(U);        
        delete [] pi;
        
}



int main(int argc , char* argv[]){
 n = stoi(argv[1]);
 num_thread = stoi(argv[2]);
 testcase = argv[3];
 if(argc > 4 && stoi(argv[4])==1){randomize=true;}
 cout<<"Size of the matrix : "<<n<<"*"<<n<<endl;
 cout<<"Total number of threads for parallel execution : "<<num_thread<<endl;
 srand48((unsigned int) time(nullptr));
 lu_decomposition();
 return 0;
}
