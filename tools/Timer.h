#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>

// TODO: There are some weird Linux specific compiler errors that arise for this struct. The solution is known although not yet ported.

struct Timer {
    Timer(const char* message, unsigned int n_tasks = 1, const char* task_name = "Task") {
        this->message = message;
        this->n_tasks = n_tasks;
        this->task_name = task_name;
        start = std::chrono::high_resolution_clock::now();
    }

    ~Timer() {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        float ms;
        if (n_tasks == 1) {
            ms = duration.count() * 1000.0f;
            std::cout << message << " " << ms << " ms\n" << std::endl;
        } else {
            ms = duration.count() * 1000.0f;
            std::cout << message << " " << ms << " ms" << std::endl;
            std::cout << task_name << " count: " << n_tasks << std::endl;
            std::cout << task_name << " average: " << ms / (float)n_tasks << " ms\n" << std::endl;
        }
    }
    const char* message;
    const char* task_name;
    unsigned int n_tasks;
    std::chrono::high_resolution_clock::time_point start, end;
    std::chrono::duration<float> duration;
};

#endif
