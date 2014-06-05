#ifndef COSMO_UTILS_TIMER_H
#define COSMO_UTILS_TIMER_H

#include "../cosmo.h"

namespace cosmo
{

/** Timer -- an individual timer **/
class Timer
{
public:
  Timer() : m_secs(0) {}

  inline double time() const { return m_secs; }

  inline void start() { clock_gettime(CLOCK_MONOTONIC, &m_starttime); }
  void stop();
  void reset();

  Timer operator+(const Timer &t2);
  Timer operator-(const Timer &t2);

private:
  double m_secs;
  struct timespec m_starttime;
  struct timespec m_stoptime;
};

/** TimerManager -- access timers via TM["my_timer"].start() **/
class TimerManager
{
public:
  TimerManager() {};

  inline Timer& operator[](std::string key) { return m_timers[key]; }

  friend std::ostream& operator<<(std::ostream &ostr, TimerManager T);

private:
  std::map<std::string, Timer> m_timers;
};

extern TimerManager _timer;

}

#endif
