
#include "Timer.h"

#include <iomanip>

namespace cosmo
{

Timer Timer::operator+(const Timer &t2)
{
  Timer t;
  t.m_secs = time() + t2.time();
  return t;
}

Timer Timer::operator-(const Timer &t2)
{
  Timer t;
  t.m_secs = time() - t2.time();
  return t;
}

void Timer::stop()
{
  clock_gettime(CLOCK_MONOTONIC, &m_stoptime);
  m_secs += (double)(m_stoptime.tv_sec - m_starttime.tv_sec);
  m_secs += (m_stoptime.tv_nsec - m_starttime.tv_nsec)*1e-9;
}

void Timer::reset()
{
  m_secs = 0.;
  m_starttime.tv_sec  = 0;
  m_starttime.tv_nsec = 0;
  m_stoptime.tv_sec   = 0;
  m_stoptime.tv_nsec  = 0;
}

std::ostream& operator<<(std::ostream &ostr, const Timer &t)
{
  ostr << std::fixed << std::setprecision(3) << t.time() << "s";
  return ostr;
}

std::ostream& operator<<(std::ostream &ostr, TimerManager T)
{
  std::map<std::string,Timer>::iterator it;

  ostr << "==== TimerManager ====" << std::endl;
  for(it = T.m_timers.begin(); it != T.m_timers.end(); ++it) {
    ostr << "  " << it->first << ": " << it->second << std::endl;
  }

  return ostr;
}

}
