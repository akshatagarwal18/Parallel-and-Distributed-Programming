#include<boost/config.hpp>
#include "mapreduce.hpp"
#include<iostream>
#include<fstream>
#include<chrono>
#include<unordered_map>

using namespace std;
using namespace std::chrono;

vector<double> page_rank;
unordered_map<int,vector<int>>adj_list;
int n; // number of nodes

namespace pr_calc{
    template<typename MapTask>
    class data_source : mapreduce::detail::noncopyable{
        private:
        int sequence_;
        public:
        data_source() : sequence_(0){}
        bool const setup_key(typename MapTask::key_type &key){
            key = sequence_++;
            return key<n;
        }
        bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value){
            value = adj_list[key];
            return true;
        }
    };

    struct map_task : public mapreduce::map_task<int,vector<int>>{
        template<typename Runtime>
        void operator() (Runtime &runtime, key_type const &key, value_type const& value) const{
            for(int i = 0 ; i < value.size() ; i++){
                runtime.emit_intermediate(value[i],page_rank[key]/adj_list[key].size());
            }
        }
    };

    struct reduce_task : public mapreduce::reduce_task<int,double>{
        template<typename Runtime, typename It>
        void operator() (Runtime& runtime, key_type const &key, It it, It ite) const{
            reduce_task::value_type result = 0;
            for(; it!=ite; it++){
                result+=(*it);
                runtime.emit(key,result);
            }
        }
    };

    typedef mapreduce::job<pr_calc::map_task, pr_calc::reduce_task,
     mapreduce::null_combiner, pr_calc::data_source<pr_calc::map_task>>job;
}

double diff(vector<double>& page_rank, vector<double>& prev){
    double d = 0.0;
    for(int i = 0 ; i < n ; i++){
        d+=(page_rank[i]-prev[i])*(page_rank[i]-prev[i]);
    } return d;
}




int main(int argc, char* argv[]){
ifstream input(argv[1]);
ofstream output;
output.open(argv[3]);
int from_page, to_page; n = 0 ;
double alpha = 0.85 , dangling, average, sum = 0.0 ;


while(!input.eof()){
    input>>from_page>>to_page;
    n = max(n,1+max(to_page,from_page));
    adj_list[from_page].push_back(to_page);
}

vector<double> pr_init(n,1.0/n);
vector<double> pr_empty(n,0.0);
vector<double> prev(n,0.0);

int iteration = 0 ;

page_rank = pr_init;
mapreduce::specification spec;
mapreduce::results result;

auto start = high_resolution_clock::now();
while(diff(page_rank,prev) > 1e-12){
    dangling = average = 0.0;
    for(int i = 0 ; i < n ; i++ ){
        average+=page_rank[i]/n;
        if(adj_list[i].size()==0){
            dangling+=page_rank[i]/n;
        }
    }

    pr_calc::job::datasource_type datasource;
    pr_calc::job mr(datasource,spec);
    mr.run<mapreduce::schedule_policy::cpu_parallel<pr_calc::job>>(result);
    prev = page_rank;
    page_rank = pr_empty;
    for(auto it = mr.begin_results(); it!= mr.end_results(); it++){
        page_rank[it->first] = it->second;
    }
    for(int i = 0 ; i < n ; i++){
        page_rank[i] = alpha*page_rank[i] + alpha*dangling + (1-alpha)*average;        
    }
}
auto end = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(end-start);

for(int i = 0 ; i < n ; i++){
    output<<i<<" = "<<page_rank[i]<<endl;
    sum+=page_rank[i];
}
cout<<duration.count()<<endl;

cout<<sum<<endl;
return 0;
}
