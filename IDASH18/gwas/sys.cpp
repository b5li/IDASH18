////////////////////////////////////////////////////////////////////////////////
// System related utilities
////////////////////////////////////////////////////////////////////////////////

#include "sys.h"

MemoryUsage getMemoryUsage() {
   FILE* file = fopen("/proc/self/status", "r");
   MemoryUsage mem = {0,0};
   char line[128];

   while(fgets(line, 128, file) != NULL) {
      if(strncmp(line, "VmRSS:", 6) == 0) {
         std::istringstream iss(line+6);
         iss >> mem.vmrss;
      } else if(strncmp(line, "VmPeak:", 7) == 0) {
         std::istringstream iss(line+7);
         iss >> mem.vmpeak;
      }
   }
   fclose(file);
   return mem;
}
