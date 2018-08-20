#ifndef SYS_H_
#define SYS_H_

#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <sys/sysinfo.h>

struct MemoryUsage {
   long vmpeak;
   long vmrss;
};

MemoryUsage getMemoryUsage();


#endif
