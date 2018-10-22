/**************************************************************************

   Copyright © 2018 Sergey Poluyan <svpoluyan@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

**************************************************************************/

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
