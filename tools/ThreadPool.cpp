//
// Created by Zach Miller on 1/14/24.
//

#include "ThreadPool.h"

#include <iostream>

void ThreadPool::Start(uint32_t n_threads, const std::function<void()>& start_task) {
    const uint32_t n_max_threads = std::thread::hardware_concurrency();  // Max # of threads the system supports
    if (n_max_threads < n_threads) {
        std::cout << "ThreadPool warning: Specified number of threads (" << n_threads << ") exceeds system maximum (" << n_max_threads << "). Setting number of threads to system maximum."
                  << std::endl;
    }
    if (n_threads <= n_max_threads) {
        for (uint32_t i = 0; i < n_threads; i++) {
            threads.emplace_back(std::thread(&ThreadPool::StartProcess, this, start_task));
        }
    } else {
        for (uint32_t i = 0; i < n_max_threads; i++) {
            threads.emplace_back(std::thread(&ThreadPool::StartProcess, this, start_task));
        }
    }
}

void ThreadPool::StartProcess(const std::function<void()>& start_task) {
    {
        std::unique_lock<std::mutex> lock(thread_state_mutex);
        thread_busy[std::this_thread::get_id()] = true;
    }
    start_task();
    {
        std::unique_lock<std::mutex> lock(thread_state_mutex);
        thread_busy[std::this_thread::get_id()] = false;
    }
    ThreadLoop();
}

void ThreadPool::ThreadLoop() {
    while (true) {
        std::function<void()> job;
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            mutex_condition.wait(lock, [this] { return !jobs.empty() || should_terminate; });
            if (should_terminate) {
                return;
            }
            job = jobs.front();
            jobs.pop();
        }
        {
            std::unique_lock<std::mutex> lock(thread_state_mutex);
            thread_busy[std::this_thread::get_id()] = true;
        }
        job();
        {
            std::unique_lock<std::mutex> lock(thread_state_mutex);
            thread_busy[std::this_thread::get_id()] = false;
            if (this->busy_cv != nullptr) {
                busy_cv->notify_one();
            }
        }
    }
}

void ThreadPool::QueueJob(const std::function<void()>& job) {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        jobs.push(job);
    }
    mutex_condition.notify_one();
}

void ThreadPool::QueueJob(const std::function<void()>& job, std::condition_variable* busy_cv) {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        jobs.push(job);
    }
    {
        std::unique_lock<std::mutex> lock(thread_state_mutex);
        this->busy_cv = busy_cv;
    }
    mutex_condition.notify_one();
}

bool ThreadPool::Busy() {
    bool poolbusy;
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        poolbusy = !jobs.empty();
        if (!poolbusy) {
            for (auto& it : thread_busy) {
                if (it.second) {
                    poolbusy = true;
                }
            }
        }
    }
    return poolbusy;
}

size_t ThreadPool::QueueSize() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    return jobs.size();
}

void ThreadPool::Stop() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        should_terminate = true;
    }
    mutex_condition.notify_all();
    for (std::thread& active_thread : threads) {
        active_thread.join();
    }
    threads.clear();
}
