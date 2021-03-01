#ifndef HPP_CPPSCRIPT
#define HPP_CPPSCRIPT

#include <thread>
#include <chrono>
#include <iostream>
namespace cppscript
{

    class _base_timer
    {
    protected:
        _base_timer() : state(new bool(true)) {}
        _base_timer(_base_timer && to)
            : mine(std::move(to.mine)), state(to.state) { to.state = nullptr; }
        _base_timer& operator= (_base_timer && to)
            {
                std::cerr << "check1a " << to.mine.get_id() << "\n";
                mine.swap(to.mine);
                std::cerr << "check1b " << mine.get_id() << "\n";
                state = to.state;
                std::cerr << "check1c\n";
                to.state = nullptr;
                std::cerr << "check1d\n";
                return *this;
            }
        ~_base_timer() { if (state) delete state; }
        bool active() const { return *state; }
        void cancel() { *state = false; }
        void watch(std::thread && t)
            { mine = std::move(t); }
        std::thread mine;
        bool* state;
    };

    class timeout : protected _base_timer
    {
        template <typename F>
        friend timeout set_timeout ( F, int );
        friend void clear_timeout ( timeout& );
    };

    class interval : protected _base_timer
    {
        template <typename F>
        friend interval set_interval ( F, int );
        friend void clear_interval ( interval& );
    };


template <typename F>
timeout set_timeout ( F func, int millis )
{
    timeout t;
    bool* active = t.state;
    t.watch( std::thread([func,millis,active](){
        std::this_thread::sleep_for( std::chrono::milliseconds(millis) );
        if (*active) { func(); }
    }) );
    return std::move(t);
}

template <typename F>
interval set_interval ( F func, int millis )
{
    interval t;
    bool* active = t.state;
    t.watch( std::thread([func,millis,active](){
        std::this_thread::sleep_for( std::chrono::milliseconds(millis) );
        while (*active) {
            func();
            std::this_thread::sleep_for( std::chrono::milliseconds(millis) );
        }
    }) );
    return std::move(t);
}

void clear_interval ( interval& t ) { t.cancel(); }
void clear_timeout ( timeout& t ) { t.cancel(); }

}

#endif
