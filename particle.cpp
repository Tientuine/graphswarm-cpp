#include "particle.hpp"

int main ()
{
    particle p1;
    particle p2 { 1, 2, 3, 4, 5 };

    p2[0] = 5.0;
    //p2.at(1) = 2.5;

    std::shared_ptr<swarmer> s1 {new swarmer};
    std::shared_ptr<swarmer> s2 {new swarmer(*s1)};

    //s1.start_swarming();
    //s1.join();
    s1->run();
    s2->run();
    s1->join();
    s2->join();

    return 0;
}
