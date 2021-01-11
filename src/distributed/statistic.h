#include <chrono>
#include <fstream>
#include "unistd.h"
#include <sys/syscall.h>
#include <algorithm>
#include <functional>

#define gettid() syscall(SYS_gettid)

using namespace std::chrono;

struct FunctionRunTime {
    std::string name;
    high_resolution_clock::time_point st;
    float time;
    int times;

    FunctionRunTime(std::string s, high_resolution_clock::time_point st, float time, int times) 
        : name(s), st(st), time(time), times(times) {}
        
    bool operator <(const FunctionRunTime& func) const {
         return time < func.time;
    }
    
    bool operator > (const FunctionRunTime& func) const {
         return time > func.time;
    }
};


struct TimeStatistics {
    
    std::vector<FunctionRunTime> func;
    high_resolution_clock::time_point app_st;
    std::mutex  mtx;

    TimeStatistics(){
        app_st = std::chrono::high_resolution_clock::now(); 
    }

    void start(const std::string &name) {
        std::string fname = name;// + " " + std::to_string(gettid());
        std::unique_lock <std::mutex> lock(mtx); 
        int size = func.size();
        int i = 0;
        for(i; i < size; i++) {
            if(func[i].name == fname) {
                func[i].st = std::chrono::high_resolution_clock::now(); 
                return;
            } 
        }
        func.emplace_back(FunctionRunTime(fname, std::chrono::high_resolution_clock::now(), 0.0, 0));
    }

    void end(const std::string &name) {
        std::string fname = name;// + " " + std::to_string(gettid());
        std::unique_lock <std::mutex> lock(mtx); 
        int size = func.size();
        int i = 0;
        for(i; i < size; i++) {
            if(func[i].name == fname) {
                func[i].time += duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - func[i].st).count();
                func[i].times++;
                return;
            } 
        }
        warn("cant find func fname ", fname, "\n");
    }
    
    void print(int frame, int rank) {
        std::ofstream os = std::ofstream("out/time_"+ std::to_string(frame) + "_" + std::to_string(rank));
        float all_time = duration_cast<std::chrono::milliseconds>(high_resolution_clock::now() - app_st).count();
        sort(func.begin(), func.end(), std::greater<FunctionRunTime>());
        os<<"| Statistic : numbers of call, time \n";
        for(int i = 0; i < func.size(); i++) {
            os<<"| "<<func[i].name<<"  ( "<<func[i].times<<" , "<<func[i].time<<" , "<<func[i].time / all_time * 100<<"% )\n";
        }     
        func.clear();
    } 

};

struct TimeStatistics statistics;


