
#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <thread>
#include <condition_variable>
#include <mutex>
#include <vector>
#include <queue>
#include <memory>
#include <exception>
#include <future>
#include <functional>

class ThreadPool {
public:
   ThreadPool(size_t n);
   ~ThreadPool();

private:
   class Task {
   private:
      ThreadPool const* pool_;
      int id_;
   public:
      Task(ThreadPool const* pool, const int id)
         : pool_(pool), id_(id) {
         // nothing to do
      }

      void run();
   };
   
   std::vector<std::thread> worker_;

   std::mutex pool_guard_;
   std::conditional_variable pool_cv_;

   bool is_active_;
};

#endif



