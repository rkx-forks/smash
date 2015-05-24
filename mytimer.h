#ifndef MYTIMER_H
#define MYTIMER_H

#include <boost/timer/timer.hpp>

#define OLD_TIMER(msg) do {\
      boost::timer::auto_cpu_timer t(msg);

#define OLD_ENDTIMER } while (0)

#define TIMER(msg) do {\
  clock_t start_t = clock();\
  string msg_s = msg;

#define ENDTIMER clock_t end_t = clock();\
  cout << msg_s << " ";\
  cout << ((float)(end_t - start_t))/CLOCKS_PER_SEC << "s" << endl;\
} while (0)

#endif
