#include <chrono>
#include <string>

#ifndef TIMER_H
#define TIMER_H

class Timer {
    typedef std::chrono::nanoseconds nanoseconds;
    typedef std::chrono::microseconds microseconds;
    typedef std::chrono::milliseconds milliseconds;
    typedef std::chrono::seconds seconds;
    typedef std::chrono::minutes minutes;
    typedef std::chrono::hours hours;
    typedef std::chrono::high_resolution_clock::time_point time_point;

    public:
        void start() {
            _start_time = std::chrono::high_resolution_clock::now();
        }

        void restart() {
            start();
        }

        std::string time_elapsed() {
            if (seconds_elapsed() < 3) {
                return std::to_string(milliseconds_elapsed()) + " milliseconds";
            }
            if (seconds_elapsed() <= 300) {
                return std::to_string(seconds_elapsed()) + " seconds";
            }
            if (minutes_elapsed() <= 300) {
                return std::to_string(minutes_elapsed()) + " minutes " +
                       std::to_string(seconds_elapsed() - 60*minutes_elapsed()) + " seconds";
            }
            if (hours_elapsed() <= 72) {
                return std::to_string(hours_elapsed()) + " hours " + 
                    std::to_string(minutes_elapsed() - 60*hours_elapsed()) + " seconds";
            }

            return std::to_string(hours_elapsed()) + " hours";
        }

        unsigned int nanoseconds_elapsed() {
            return std::chrono::duration_cast<nanoseconds>(_elapsed()).count();
        }

        unsigned int microseconds_elapsed() {
            return std::chrono::duration_cast<microseconds>(_elapsed()).count();
        }

        unsigned int milliseconds_elapsed() {
            return std::chrono::duration_cast<milliseconds>(_elapsed()).count();
        }

        unsigned int seconds_elapsed() {
            return std::chrono::duration_cast<seconds>(_elapsed()).count();
        }

        unsigned int minutes_elapsed() {
            return std::chrono::duration_cast<minutes>(_elapsed()).count();
        }

        unsigned int hours_elapsed() {
            return std::chrono::duration_cast<hours>(_elapsed()).count();
        }

    private:
        time_point _start_time;
        nanoseconds _elapsed() {
            return std::chrono::high_resolution_clock::now() - _start_time;
        }
};
#endif
