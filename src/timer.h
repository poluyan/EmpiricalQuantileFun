#ifndef TIME_H
#define TIME_H

#include <ctime>
#include <chrono>

namespace timer
{
    class Timer
    {
    public:
        Timer(): beg_(clock_::now()) {}
        void reset();
        double elapsed_nanoseconds() const;
        double elapsed_microseconds() const;
        double elapsed_milliseconds() const;
        double elapsed_seconds() const;
        double elapsed_minutes() const;
        double elapsed_hours() const;
    private:
        typedef std::chrono::high_resolution_clock clock_;
        typedef std::chrono::duration<double, std::nano> nanoseconds_;
        typedef std::chrono::duration<double, std::micro> microseconds_;
        typedef std::chrono::duration<double, std::milli> milliseconds_;
        typedef std::chrono::duration<double, std::ratio<1> > seconds_;
        typedef std::chrono::duration<double, std::ratio<60> > minutes_;
        typedef std::chrono::duration<double, std::ratio<3600> > hours_;
        std::chrono::time_point<clock_> beg_;
    };
}

#endif
