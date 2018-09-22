#include "timer.h"

namespace timer
{
    void Timer::reset()
    {
        beg_ = clock_::now();
    }
    double Timer::elapsed_nanoseconds() const
    {
        return std::chrono::duration_cast<nanoseconds_>(clock_::now() - beg_).count();
    }
    double Timer::elapsed_microseconds() const
    {
        return std::chrono::duration_cast<microseconds_>(clock_::now() - beg_).count();
    }
    double Timer::elapsed_milliseconds() const
    {
        return std::chrono::duration_cast<milliseconds_>(clock_::now() - beg_).count();
    }
    double Timer::elapsed_seconds() const
    {
        return std::chrono::duration_cast<seconds_>(clock_::now() - beg_).count();
    }
    double Timer::elapsed_minutes() const
    {
        return std::chrono::duration_cast<minutes_>(clock_::now() - beg_).count();
    }
    double Timer::elapsed_hours() const
    {
        return std::chrono::duration_cast<hours_>(clock_::now() - beg_).count();
    }
}
