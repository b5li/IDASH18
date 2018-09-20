#ifndef SYS_H_
#define SYS_H_

#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <chrono>
#include <sys/types.h>
#include <sys/sysinfo.h>
struct MemoryUsage {
   long vmpeak;
   long vmrss;
};

MemoryUsage getMemoryUsage();

class TimerAnchor {
   char const * title;
   std::chrono::steady_clock::time_point start;
public:

   TimerAnchor(const char * _title) : title(_title), start(std::chrono::steady_clock::now()) { }

   void printTimeInfo() const;

   ~TimerAnchor() {
      printTimeInfo();
   }
};



#endif
