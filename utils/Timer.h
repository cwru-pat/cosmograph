#ifndef COSMO_UTILS_TIMER_H
#define COSMO_UTILS_TIMER_H

/** Timer -- an individual timer **/
class Timer
{
public:
  Timer() : m_secs(0) {}

  inline Timer operator+(const Timer &t2)
  {
    Timer t;
    t.m_secs = time() + t2.time();
    return t;
  }

  inline double time() const
  {
    return m_secs;
  }

  inline void start()
  {
    clock_gettime(CLOCK_MONOTONIC, &m_starttime);
  }

  inline void stop()
  {
    clock_gettime(CLOCK_MONOTONIC, &m_stoptime);
    m_secs += (double)(m_stoptime.tv_sec - m_starttime.tv_sec);
    m_secs += (m_stoptime.tv_nsec - m_starttime.tv_nsec)*1e-9;
  }

  inline void reset()
  {
    m_secs = 0.;
    m_starttime.tv_sec  = 0;
    m_starttime.tv_nsec = 0;
    m_stoptime.tv_sec   = 0;
    m_stoptime.tv_nsec  = 0;
  }

  double m_secs;
private:
  struct timespec m_starttime;
  struct timespec m_stoptime;
};

inline std::ostream& operator<<(std::ostream &ostr, const Timer &t)
{
  ostr << std::fixed << std::setprecision(3) << t.time() << "s";
  return ostr;
}

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

std::ostream& operator<<(std::ostream &ostr, TimerManager T)
{
  std::map<std::string,Timer>::iterator it;

  ostr << "==== TimerManager ====" << std::endl;
  for(it = T.m_timers.begin(); it != T.m_timers.end(); ++it) {
    ostr << "  " << it->first << ": " << it->second << std::endl;
  }

  return ostr;
}

#endif
