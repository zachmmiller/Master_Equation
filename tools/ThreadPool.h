//
// Created by Zach Miller on 1/14/24.
//

#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <complex>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <unordered_map>
#include <vector>

class ThreadPool {
   public:
    void Start(uint32_t n_threads, const std::function<void()>& start_task);
    void QueueJob(const std::function<void()>& job);
    void QueueJob(const std::function<void()>& job, std::condition_variable* busy_cv);
    void Stop();
    bool Busy();
    size_t QueueSize();

   private:
    void StartProcess(const std::function<void()>& start_task);
    void ThreadLoop();

    bool should_terminate = false;            // Tells threads to stop looking for jobs
    std::mutex queue_mutex;                   // Prevents data races to the job queue
    std::condition_variable mutex_condition;  // Allows threads to wait on new jobs or termination
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> jobs;

    std::mutex thread_state_mutex;
    std::unordered_map<std::thread::id, bool> thread_busy;
    std::condition_variable* busy_cv;
};

#endif  // THREADPOOL_H
