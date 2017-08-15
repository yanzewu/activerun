#pragma once
#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>


class FixedThreadPool {
public:

	FixedThreadPool() {

	}

	void init(size_t capacity, bool auto_finish=true) {
		this->auto_finish = auto_finish;
		interrputed = false;
		execute_count = 0;
		pool.reserve(capacity);
		for (size_t i = 0; i < capacity; i++) {
			pool.emplace_back(std::bind(&FixedThreadPool::start_job, this, i));
		}
	}

	~FixedThreadPool() {
		{
			std::unique_lock<std::mutex> l(lock_start);
			interrputed = true;
			condvar_start.notify_all();
		}
		for (auto& t : pool) {
			t.join();
		}
	}

	void interrupt() {
		interrputed = true;
		condvar_start.notify_all();
	}

	void submit(std::function<void(void)> func, bool start_rightnow=true) {

		{
			std::unique_lock<std::mutex> l(lock_finish);
			execute_count++;
			condvar_finish.notify_one();
		}
		{
			std::unique_lock<std::mutex> l(lock_start);
			job_queue.emplace(std::move(func));
			if (start_rightnow)
				condvar_start.notify_one();
		}
	}

	void wait() {
		std::unique_lock<std::mutex> l(lock_finish);
		while (!interrputed && execute_count > 0) {
			condvar_finish.wait(l);
		}
	}

	void start_all() {
		condvar_start.notify_all();
	}

	void start_job(size_t index) {
		std::function<void(void)> job;
		while (1) {
			// wait until acquire job
			{
				std::unique_lock<std::mutex> l(lock_start);
				while (!interrputed && (auto_finish && job_queue.empty())) {
					condvar_start.wait(l);
				}
				
				if (auto_finish && job_queue.empty()) {
					return;
				}

				job = std::move(job_queue.front());
				job_queue.pop();
			}

			job();
			
			{
				std::unique_lock<std::mutex> l(lock_finish);
				execute_count--;
				condvar_finish.notify_one();
			}
		}
	}

	std::mutex lock_start, lock_finish;
	std::condition_variable condvar_start, condvar_finish;
	std::vector<std::thread> pool;
	std::queue<std::function<void(void)> > job_queue;
	int execute_count;
	bool interrputed;
	bool auto_finish;
};

#endif // !THREADPOOL_H
