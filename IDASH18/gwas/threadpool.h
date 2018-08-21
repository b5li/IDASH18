
#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <thread>
#include <condition_variable>
#include <mutex>
#include <vector>
#include <memory>
#include <exception>
#include <assert.h>
#include <future>
#include <functional>

#define DEFAULT_NUM_THREADS 4

#define LogicError(s)                           \
   assert(0&&(s))

#define ResourceError(s)                        \
   assert(0&&(s))

struct PartitionInfo {
   long nintervals;  // number of intervals
   long intervalsz;  // interval size
   long nsintervals; // number of small intervals

   explicit
   PartitionInfo(long sz, long nt = DEFAULT_NUM_THREADS) 
   // partitions [0..sz) into nintervals intervals,
   // so that there are nsintervals of size intervalsz-1
   // and nintervals-nsintervals of size intervalsz
   {
      if (sz <= 0) {
         nintervals = intervalsz = nsintervals = 0;
         return;
      }

      if (nt <= 0) LogicError("PartitionInfo: bad args");

      if (sz < nt) {
         nintervals = sz;
         intervalsz = 1;
         nsintervals = 0;
         return;
      }

      nintervals = nt;

      long q, r;
      q = sz/nt;
      r = sz - nt*q;

      if (r == 0) {
         intervalsz = q;
         nsintervals = 0;
      }
      else {
         intervalsz = q+1;
         nsintervals = nt - r;
      }
   }

   long NumIntervals() const { return nintervals; }

   void interval(long& first, long& last, long i) const
   // [first..last) is the ith interval -- no range checking is done
   {

#if 1
      // this is the logic, naturally expressed
      if (i < nsintervals) {
         first = i*(intervalsz-1);
         last = first + (intervalsz-1);
      }
      else {
         first = nsintervals*(intervalsz-1) + (i-nsintervals)*intervalsz;
         last = first + intervalsz;
      }
#else
      // this is the same logic, but branch-free (and portable)
      // ...probably unnecessary optimization
      
      long mask = -long(reinterpret_cast<unsigned long>(i-nsintervals) >> (NTL_BITS_PER_LONG-1));
      // mask == -1 if i < nsintervals, 0 o/w
 
      long lfirst = i*(intervalsz-1);
      lfirst += long((~reinterpret_cast<unsigned long>(mask)) & cast_unsigned(i-nsintervals));
      // lfirst += max(0, i-nsintervals)

      long llast = lfirst + intervalsz + mask;

      first = lfirst;
      last = llast;
#endif
   }

};

class ThreadPool {

// public:
   // ThreadPool(size_t n);
   // ~ThreadPool();

private:

   template<class T>
   class SimpleSignal {
   private:
     T val; 
     std::mutex m;
     std::condition_variable cv;
   
     SimpleSignal(const SimpleSignal&); // disabled
     void operator=(const SimpleSignal&); // disabled
   
   public:
     SimpleSignal() : val(0) { }
   
     T wait() 
     {
       std::unique_lock<std::mutex> lock(m);
       cv.wait(lock, [&]() { return val; } );
       T old_val = val;
       val = 0;
       return old_val;
     }
   
     void send(T new_val)
     {
       std::lock_guard<std::mutex> lock(m);
       val = new_val;
       cv.notify_one();
     }
   };
   
   
   template<class T, class T1>
   class CompositeSignal {
   private:
     T val; 
     T1 val1;
     std::mutex m;
     std::condition_variable cv;
   
     CompositeSignal(const CompositeSignal&); // disabled
     void operator=(const CompositeSignal&); // disabled
   
   public:
     CompositeSignal() : val(0) { }
   
     T wait(T1& _val1) 
     {
       std::unique_lock<std::mutex> lock(m);
       cv.wait(lock, [&]() { return val; } );
       T _val = val;
       _val1 = val1;
       val = 0;
       return _val;
     }
   
     void send(T _val, T1 _val1)
     {
       std::lock_guard<std::mutex> lock(m);
       val = _val;
       val1 = _val1;
       cv.notify_one();
     }
   };

   class Task {
   private:
      ThreadPool * pool_;
      // int id_;
   public:
      Task(ThreadPool * pool)
         : pool_(pool) {
         // nothing to do
      }
      ThreadPool * getThreadPool() { return pool_; }
      virtual void run(long index) = 0;
   };


   class TerminateTask : public Task {
   public:
      TerminateTask() : Task(0) { }
      void run(long) { }
   };

   template<class F>
   class FunctionTask : public Task {
   public:
      F const& f;
      PartitionInfo const& pinfo;
      FunctionTask(ThreadPool * pool, F const& _f, PartitionInfo const& _pinfo)
         : Task(pool), f(_f), pinfo(_pinfo) {
         // nothing to do
      }
      void run(long index) {
         long first, last;
         pinfo.interval(first, last, index);
         f(first, last);
      }
   };

   struct AutomaticThread {
      CompositeSignal< Task *, long > localSignal;
      TerminateTask term;
      std::thread t;
   
      AutomaticThread() : t(worker, &localSignal) { 
         // std::cerr << "starting thread " << t.get_id() << "\n";
      }
   
      ~AutomaticThread() {
        // std::cerr << "stopping thread " << t.get_id() << "...";
        localSignal.send(&term, -1);
        t.join();
        // std::cerr << "\n";
      }
   };

   std::vector<AutomaticThread *> thread_;

   // std::mutex pool_guard_;
   // std::condition_variable pool_cv_;

   SimpleSignal<bool> globalSignal_;

   bool is_active_;             // FixMe: should this be atomic<bool> ?

   std::atomic<long> counter_;

   ThreadPool(const ThreadPool&) = delete;
   void operator=(const ThreadPool&) = delete;

   void launch(Task *task, long index) {
      thread_[index-1]->localSignal.send(task, index);
      // we use threadVec[index-1] to allow for the fact
      // that we want the current thread to have index 0
   }

   void begin(long cnt) {
      is_active_ = true;
      counter_ = cnt;
   }

   void end() {
      globalSignal_.wait();
      is_active_ = false;
   }

   static void runOneTask(Task *task, long index) {
      ThreadPool * pool = task->getThreadPool();

      task->run(index);

      if(--(pool->counter_) == 0) {
         pool->globalSignal_.send(true);
      }
   }

   static void worker(CompositeSignal< Task *, long > * localSignal) {
      for (;;) {
         long index = -1;
         Task *task = localSignal->wait(index);
         if(index == -1) {
            return;
         }
         runOneTask(task, index);
      }
   }

public:

   long numThreads() const {
      return thread_.size();
   }
   bool active() const {
      return is_active_;
   }

   explicit ThreadPool(size_t n)
      : is_active_(false), counter_(0) {
      for(int i = 1; i < n; i++) { // create n-1 threads
         AutomaticThread * t = new AutomaticThread();
         thread_.push_back(t);
      }
   }
         
   ~ThreadPool() {
      if(active()) {
         assert(0&&"Destructed while active");
      }
      while(!thread_.empty()) {
         AutomaticThread * t = thread_.back();
         thread_.pop_back();
         delete t;
      }
   }

   template<class Fct>
   void exec_range(long sz, const Fct& fct) {
      if(active()) {
         LogicError("ThreadPool: illegal operation while active");
      }
      if(sz <= 0) {
         return;
      }

      PartitionInfo pinfo(sz, numThreads());
      long cnt = pinfo.NumIntervals();

      FunctionTask<Fct> task(this, fct, pinfo);

      begin(cnt);
      for(long t = 1; t < cnt; t++) {
         launch(&task, t);
      }
      runOneTask(&task, 0);
      end();
   }
};

template<class Fct>
static void relaxed_exec_range(ThreadPool *pool, long sz, const Fct& fct) {
   if(sz <= 0) {
      return;
   }
   if(!pool || pool->active() || sz == 1) {
      fct(0, sz);
   } else {
      pool->exec_range(sz, fct);
   }
}


extern
ThreadPool *threadPool_ptr__;

inline ThreadPool *getThreadPool() {
   return threadPool_ptr__;
}

void initThreadPool(size_t n);


inline long availableThreads() {
   ThreadPool *pool = getThreadPool();
   if(!pool || pool->active()) {
      return 1;
   } else {
      return pool->numThreads();
   }
}

#define TP_EXEC_RANGE(n, first, last)                                  \
   {                                                                   \
   relaxed_exec_range(getThreadPool(), (n),                            \
      [&](long first, long last) {                                     \


#define TP_EXEC_RANGE_END                          \
      } );                                         \
   }                                               \


#endif



