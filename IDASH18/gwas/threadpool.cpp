#include "threadpool.h"

namespace IDASH {

ThreadPool *threadPool_ptr__;

void initThreadPool(size_t n) {
   threadPool_ptr__ = new ThreadPool(n);
}

}
