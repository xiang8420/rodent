#include <chrono>
#include <fstream>
#include "unistd.h"
#include <sys/syscall.h>
#define gettid() syscall(SYS_gettid)

using namespace std::chrono;
struct TimeStatistics {

    std::vector<std::string>                       func_name;
    std::vector<high_resolution_clock::time_point> func_st;
    std::vector<float>                             func_time;
    std::vector<int>                               func_num_run;
    high_resolution_clock::time_point              app_st;

    std::mutex  mtx;

    TimeStatistics(){
        app_st = std::chrono::high_resolution_clock::now(); 
    }

    void start(const std::string &name) {
        std::string fname = name + " " + std::to_string(gettid());
        std::unique_lock <std::mutex> lock(mtx); 
        std::vector<std::string>::iterator iter = find(func_name.begin(), func_name.end(), fname);
        if( iter != func_name.end() ) {
            int pos = distance(func_name.begin(), iter);
            func_st[pos] = std::chrono::high_resolution_clock::now(); 
            func_num_run[pos]++;
        } else {
            func_name.emplace_back(fname); 
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
        std::string fname = name + " " + std::to_string(gettid());
        std::unique_lock <std::mutex> lock(mtx); 
        std::vector<std::string>::iterator iter = find(func_name.begin(), func_name.end(), fname);
        for(int i = 0; i < func_name.size(); i++) {
            if(func_name[i] == fname) {
                func_time[i] += duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - func_st[i]).count();
                return ;
            }
        }
        warn("cant find func fname\n");
        warn(fname);
    }
    
    void print(std::ofstream &os) {
        float all_time = duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - app_st).count();
        os<<"| Statistic : numbers of call, time \n";
        for(int i = 0; i < func_name.size(); i++) {
            os<<"| "<<func_name[i]<<"  ( "<<func_num_run[i]<<" , "<<func_time[i] / all_time * 100<<"% )\n";
        }     
    
    } 

};

struct TimeStatistics statistics;


