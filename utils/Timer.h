#ifndef COSMO_UTILS_TIMER_H
#define COSMO_UTILS_TIMER_H

#include <ctime>
#include <iostream>
#include <map>

namespace cosmo
{

/**
 * @brief Individual timer classes used by the TimerManager class.
 */
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

/**
 * @brief TimerManager class containing multiple timers;
 * access individual timers via, eg, TM["my_timer"].start()
 */
class TimerManager
{
public:
  TimerManager() {};

  std::string getStateString();

  inline Timer& operator[](std::string key) { return m_timers[key]; }

  friend std::ostream& operator<<(std::ostream &ostr, TimerManager T);

private:
  std::map<std::string, Timer> m_timers;
};

}

#endif
