#include "objs.hpp"

int main (int argc, char** argv)
{
    int n = atoi(argv[1]);

    std::forward_list<double> x (n,2.5);

    rosenbrock r;
    double c = r(x.begin(), x.end());
    std::cout << std::setprecision(12) << c << std::endl;
}
