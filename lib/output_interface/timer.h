#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <ctime>
#include <string>

class Timer {
  typedef std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> time_point;

public:
  Timer() { start_timer(); }

  static inline std::string get_time() {

    std::time_t result = std::time(nullptr);
    const std::tm *t = std::localtime(&result);
    return "[" + std::to_string(t->tm_hour) +
           ":" + std::to_string(t->tm_min) +
           ":" + std::to_string(t->tm_sec) + "]: ";
  }

  inline void start_timer() {
    m_time = get_cur_time();
  }

  inline double get_delta_time_sec() {
    return static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(get_cur_time() - m_time).count()) / 1000.;
  }

  inline int64_t get_delta_time_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(get_cur_time() - m_time).count();
  }

  inline int64_t get_delta_time_ns() {
    return (get_cur_time() - m_time).count();
  }

private:
  inline time_point get_cur_time() {
    return std::chrono::steady_clock::now();
  }

private:
  time_point m_time;
};

#endif //! TIMER_H