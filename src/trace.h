/*==============================================================================*
 * TRACE
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-08-24
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <omp.h>

#define TRACE_ME TraceCaller _TRACE_OBJECT(__func__);
#define TRACE(x) TraceCaller _TRACE_OBJECT(x);

double _walltime() {
    return omp_get_wtime();
}

struct TraceElement {
    std::string name;
    int num_calls = 0;
    double total_time = 0;
    double child_time = 0;
    double start_time = 0;
};

class Trace {
    private:
        std::string name;
        std::vector<TraceElement*> callstack;
        std::unordered_map<std::string, TraceElement> function_list;

    public:
        Trace(std::string name_): name(std::move(name_)) {

        };

        void write_callstack() {
            std::cout << "Callstack:" << std::endl;
            for(const auto &call: callstack) {
                std::cout << "    " << call->name << "\t" << call->total_time << std::endl;
            }
        };

        void write_profile(std::string filename) {
            std::ofstream trace_file;
            trace_file.open(filename);
            trace_file << "function_name"<< ","
                       << "total_time" << ","
                       << "exclusive_time" << ","
                       << "number_of_calls"
                       << std::endl;

            for(const auto& func: function_list) {
                if (!func.second.start_time) {
                    trace_file << func.first << ","
                               << func.second.total_time << ","
                               << func.second.total_time - func.second.child_time << ","
                               << func.second.num_calls
                               << std::endl;
                }
            }
            trace_file.close();
        };

        void enter(std::string func_name) {
            // if (!function_list.contains(func_name)) {
            if (function_list.find(func_name) == function_list.end()) {
                function_list.insert(std::make_pair(func_name, TraceElement()));
                function_list[func_name].name = func_name;
            }

            auto start_time = _walltime();
            callstack.push_back(&function_list[func_name]);
            function_list[func_name].num_calls += 1;
            function_list[func_name].start_time = start_time;
        };

        void exit(std::string func_name) {
            auto end_time = _walltime();
            auto total_time = end_time - function_list[func_name].start_time;
            function_list[func_name].total_time += total_time;
            function_list[func_name].start_time = 0;
            callstack.pop_back();
            if (!callstack.empty()) callstack.back()->child_time += total_time;
        }

};

namespace trace {
    static const int enabled = 1;
    Trace current("__TRACE_BASE__");
}

struct TraceCaller {
    std::string name;
    TraceCaller(std::string name_) : name(std::move(name_)) { if (trace::enabled) trace::current.enter(name) ;};
    ~TraceCaller() { if (trace::enabled) trace::current.exit(name); };

    void add(std::string name_);
};

#endif /* !MAIN_H */
