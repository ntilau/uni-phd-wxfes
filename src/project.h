#ifndef PROJECT_H
#define PROJECT_H

#include <string>
#include <time.h>

#include "model.h"

class timer {
public:
    timer() {
        tic();
    }
    ~timer() {
    }
    void tic() {
        lc = clock();
    }
    double toc() {
        clock_t cc = clock();
        return (cc - lc) / 1000.0;
    }
    std::string strtoc() {
        clock_t cc = clock();
        std::stringstream timing;
        timing << (cc - lc) / 1000.0 << " s";
        return timing.str();
    }
private:
    clock_t lc;
};

class project {
public:
    enum task_to_perform {
        NONE = 0,
        LOAD_FES,
        SAVE_FES,
        LOAD_HFSS,
        LOAD_AEDT,
        LOAD_POLY,
        SAVE_POLY,
        LOAD_STL,
        RUN_TETGEN,
        RUN_TRIANGLE,
        REFINE_HOMOGENEOUSLY,
        ANALYZE
    } task;
    project();
    ~project();
    mdl_core model; // all the model information is stored here
    std::string name = "";
    std::string data_path = "./";
    std::string aux_path = "./";
    std::string full_path_name = "./";
    std::string get_stats(timer&); // Statistics
    void execute_task();
    std::string get_info(); // Computer name, CPU cores and OMP threads number
    int get_num_proc(); // Number of processors
    std::vector<std::string> get_mac();
    bool check_mac();
    std::string set_priority(unsigned int);
    std::string get_priority();
    std::string get_proc_mem();
    std::string get_sys_mem();
    std::string get_loc_time();
    std::string get_gmt();
private:
    std::string get_var(const std::string name);
    int get_int(const std::string name);
};

#endif // PROJECT_H
