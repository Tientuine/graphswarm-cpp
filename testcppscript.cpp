#include "cppscript.hpp"

using namespace cppscript;

#include <iostream>
#include <iomanip>

void sayHi ()
{
    std::cout << "Hi!\n";
}

void sayHello ()
{
    std::cout << "Hello!\n";
}

void sayBye ()
{
    std::cout << "Bye!\n";
}


int main()
{
    timeout t1 = set_timeout(sayHi, 500);
    interval t2 = set_interval(sayBye, 800);
    timeout t3 = set_timeout(sayHello, 2000);

    std::this_thread::sleep_for(std::chrono::seconds(1));

    clear_timeout(t3);

    std::this_thread::sleep_for(std::chrono::seconds(4));

    std::cerr << "check1\n";
    t3 = set_timeout(sayHi, 200);
    std::cerr << "check2\n";

    clear_interval(t2);

    return 0;
}
