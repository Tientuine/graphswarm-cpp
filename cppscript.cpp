#ifndef HPP_CPPSCRIPT
#define HPP_CPPSCRIPT

#include <thread>
#include <chrono>

namespace cppscript
{

template <typename F>
std::thread set_timeout ( F func, int millis )
{
    return std::thread([](){
        std::sleep_for( std::chrono::milliseconds(millis) );
        if (this->joinable()) { func(); }
    });
}

template <typename F>
std::thread set_timeout ( F func, int millis )
{
    return std::thread([](){
        while (true) {
            std::sleep_for( std::chrono::milliseconds(millis) );
            if (this->joinable()) { func(); }
            else break;
        }
    });
}

void clear_interval ( std::thread t ) { t.detach(); }
void clear_timeout ( std::thread t ) { t.detach(); }

}

#endif
