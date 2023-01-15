#include "mpi.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include "sys/stat.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include<bits/stdc++.h>

using namespace std;
#define MASTER 0

struct MY_STRUCT{
    unordered_map<int,vector<int>>adj_list;
    vector<int>node_list;
    vector<double> page_rank;
    int num_processes;
    int graph_size;
    int num;
    int val;
};

void mapper(int rank, MAPREDUCE_NS::KeyValue *kv, void* mystr){
    MY_STRUCT *S = (MY_STRUCT *)mystr;
    int graph_size = S->graph_size;
    int num_processes = S->num_processes;
    int start = rank*(graph_size/num_processes);
    int end = (rank+1)*(graph_size/num_processes);
    if(rank == S->num_processes-1){end = graph.size();}
    for(int i = start ; i < end ; i++){

    }
}

int main(int argc, char* argv[]){
    ifstream input;
    input.open(argv[1]);
    ofstream output;
    output.open(argv[2]);

    int from_page, to_page;
    int rank, num_processes;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_Size(MPI_COMM_WORLD, &num_processes);
    if(num_processes < 1){
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(1);
    }

    MAPREDUCE_NS::MapReduce *map_red = new MAPREDUCE_NS::MapReduce(MPI_COMM_WORLD);
    vector<double>page_rank;
    string line;
    if(input.is_open()){
        while(getline(input,line)){
            istringstream input_line(line);
            input_line>>from_page>>to_page;
            adj_list[from_page].push_back(to_page);
;        }
    }

    int graph_size = adj_list.size();
    auto it ;
    for(it = adj_list.begin(); it!= adj_list.end(); it++){
        page_rank[it->first] = 1.0/(double)graph_size;
    }

    MY_STRUCT my_struct;
    my_struct.page_rank = page_rank;
    my_struct.adj_list = adj_list;
    my_struct.num_processes = num_processes;
    my_struct.graph_size = graph_size;
    my_struct.num = 0;
    my_struct.val = 0;
    MPI_Barries(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    for(int i = 0 ; i < 20 ; i++){
        map_red->map(num_processes, &mapper, &my_struct);
        map_red->convert();
        map_red->reduce(reducer, &my_struct);
        map_red->gather(1);
        delete mr;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double stop = MPI_Wtime();
    double time = stop - start;

    if(rank == MASTER)
    {
        for(int i=0; i<my_struct.graph_size; i++)
        {
            output << i << " = " << my_struct.page_rank[i] << endl;
        }
        printf("The time taken for %d graph size and %d processes using Blocking P2P Communication is %0.4fs\n"
                , graph_size, num_processes, time);
    }

    MPI_Finalize();


    return 0;
}