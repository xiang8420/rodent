#include <chrono>
#include <fstream>

using namespace std::chrono;
struct TimeStatistics {

    std::vector<std::string>                       func_name;
    std::vector<high_resolution_clock::time_point> func_st;
    std::vector<float>                             func_time;
    std::vector<int>                               func_num_run;

    TimeStatistics(){}

    void start(const std::string &name) {
        std::vector<std::string>::iterator iter = find(func_name.begin(), func_name.end(), name);
        if( iter != func_name.end() ) {
            int pos = distance(func_name.begin(), iter);
            func_st[pos] = std::chrono::high_resolution_clock::now(); 
            func_num_run[pos]++;
        } else {
            func_name.emplace_back(name); 
            func_st.emplace_back(std::chrono::high_resolution_clock::now());
            func_time.emplace_back(0.0);
            func_num_run.emplace_back(1); 
        }

       // for(int i = 0; i < func_name.size(); i++) {
       //     if(func_name[i] == name) {
       //         func_st[i]
       //         break;
       //     }
       // }

       // if()
    }

    void end(const std::string &name) {
        std::vector<std::string>::iterator iter = find(func_name.begin(), func_name.end(), name);
        if( iter != func_name.end() ) {
            int pos        = distance(func_name.begin(), iter);
            func_time[pos] += duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - func_st[pos]).count();
        } else {
            error(name);
        }
    }
    
    void print(std::ofstream &os) {
        os<<"| Statistic : numbers of call, time \n";
        for(int i = 0; i < func_name.size(); i++) {
            os<<"| "<<func_name[i]<<"  ( "<<func_num_run[i]<<" , "<<func_time[i]<<" )\n";
        }     
    
    } 

};

struct TimeStatistics statistics;


