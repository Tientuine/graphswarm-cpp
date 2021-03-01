//#include "objective.hpp"
#include "runnables.hpp"
#include <algorithm>
#include <forward_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <system_error>
#include <vector>

struct objective
{
    using param = std::forward_list<double>::const_iterator;
    using domain_type = std::pair<double,double>;

    virtual auto operator() ( param a, param b ) const -> double = 0;
    virtual auto domain ( unsigned i ) const -> domain_type = 0;
    virtual auto extremum ( unsigned i ) const -> double = 0;
};

class sphere : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto cost = std::inner_product(a, b, a, 0.0);
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-5.12, 5.12); }
    auto extremum ( unsigned i ) const -> double { return 0.0; }
};

class rosenbrock : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto c = a;
        auto cost = std::inner_product(++a, b, a, 0.0,
            std::plus<double>(),
            [](double x2, double x1) {
                auto t1 = x1 * x1  - x2;
                auto t2 = x1 - 1.0;
                return 100.0 * t1 * t1 + t2 * t2;
            }
        );
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-5.0, 10.0); }
    auto extremum ( unsigned i ) const -> double { return (i == 0 ? 0.0 : 1.0); }
};

class rastrigin : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto cost = std::accumulate(a, b, 0.0,
            [](double sum, double x) {
                static auto const TWOPI = 8.0 * std::atan(1.0);
                return sum + x * x - 10.0 * std::cos(TWOPI * x);
            }
        );            
        unsigned n = distance(a,b);
        return 10.0 * n + cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-5.12, 5.12); }
    auto extremum ( unsigned i ) const -> double { return 0.0; }
};

class griewangk : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto cost1 = std::accumulate(a, b, 0.0,
            [](double sum, double x) { return sum + x * x  / 4000.0; }
        );
        auto i = 0.0;
        auto cost2 = std::accumulate(a, b, 1.0,
            [&i](double prod, double x) {
                return prod * std::cos(x / std::sqrt(++i));
            }
        );
        return cost1 - cost2 + 1.0;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-600.0, 600.0); }
    auto extremum ( unsigned i ) const -> double { return 0.0; }
};

class shaffer_f6 : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto x1 = *(a++), x2 = *a;
        auto h = x1 * x1 + x2 * x2;
        auto denom = 1 + 0.001 * h;
        auto numer = std::sin(std::sqrt(h));  
        auto cost = 0.5 + (numer * numer - 0.5) / (denom * denom);
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-10.0, 10.0); } //!!
    auto extremum ( unsigned i ) const -> double { return 0.0; }
};

class shaffer_f6_inv : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto x1 = *(a++), x2 = *a;
        auto h = x1 * x1 + x2 * x2;
        auto denom = 1 + 0.001 * h;
        auto numer = std::sin(std::sqrt(h));  
        auto cost = 0.5 - (numer * numer - 0.5) / (denom * denom);
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-10.0, 10.0); } //!!
    auto extremum ( unsigned i ) const -> double { return 0.0; }
};

class beale : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        if (std::distance(a, b) != 2)
            { throw std::logic_error("must have exactly 2 dimensions"); }
        auto& x1 = *(a++);
        auto& x2 = *a;
        auto t1 = 1.5 - x1 * (1 - x2);
        auto t2 = 2.25 - x1 * (1 - x2 * x2);
        auto t3 = 2.625 - x1 * (1-x2*x2*x2);
        double cost = t1 * t1 + t2 * t2 + t3 * t3;
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-4.5, 4.5); }
    auto extremum ( unsigned i ) const -> double
        { return (i == 0 ? 0.0 : (i == 1 ? 3.0 : 0.5)); }
};

/*
class bohachevsky1
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        if (std::distance(a, b) != 2)
            { throw std::logic_error("must have exactly 2 dimensions"); }
        auto& x1 = *(a++);
        auto& x2 = *a;
        auto t1 = x1 * x1;
        auto t2 = 2 * x2 * x2;
        auto t3 = 0.3 * std::cos( THREE_PI * x1 + FOUR_PI * x2 );
        double cost = t1 * t1 + t2 * t2 + t3 * t3;
        return cost;
    }
    //auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-600.0, 600.0); }
    //auto extremum ( unsigned i ) const -> double { return 0.0; }
};
*/

class booth : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        if (std::distance(a, b) != 2)
            { throw std::logic_error("must have exactly 2 dimensions"); }
        auto& x1 = *(a++);
        auto& x2 = *a;
        auto t1 = x1 + 2 * x2 - 7;
        auto t2 = 2 * x1 + x2 - 5;
        double cost = t1 * t1 + t2 * t2;
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-10.0, 10.0); }
    auto extremum ( unsigned i ) const -> double
        { return (i == 0 ? 0.0 : (i == 1 ? 1.0 : 3.0)); }
};

class branin : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        if (std::distance(a, b) != 2)
            { throw std::logic_error("must have exactly 2 dimensions"); }
        auto& x1 = *(a++);
        auto& x2 = *a;
        auto t1 = 1.5 - x1 * (1 - x2);
        auto t2 = 2.25 - x1 * (1 - x2 * x2);
        auto t3 = 2.625 - x1 * (1-x2*x2*x2);
        double cost = t1 * t1 + t2 * t2 + t3 * t3;
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type
        { return (i == 1 ? std::make_pair(-5.0, 10.0) : std::make_pair(-0.0, 15.0)); }
    auto extremum ( unsigned i ) const -> double
        { return (i == 0 ? 0.397887 : (i == 1 ? 9.42478 : 2.475)); }
};

class colville : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        if (std::distance(a, b) != 4)
            { throw std::logic_error("must have exactly 4 dimensions"); }
        auto& x1 = *(a++);
        auto& x2 = *(a++);
        auto& x3 = *(a++);
        auto& x4 = *a;
        auto t1 = x1 * x1 - x2;
        auto t2 = x1 - 1;
        auto t3 = x3 - 1;
        auto t4 = x3 * x3 - x4;
        auto t5 = x4 - 1;
        auto t6 = x2 - 1;
        double cost = 100.0 * t1 * t1 + t2 * t2 + t3 * t3 + 90.0 * t4 * t4
            + 10.1 * (t3 * t3 + t5 * t5) + 19.8 * t6 * t5;
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-10.0, 10.0); }
    auto extremum ( unsigned i ) const -> double { return (i == 0 ? 0.0 : 1.0); }
};

class dixon_price : public objective
{
public:
    auto operator() ( param a, param b ) const -> double
    {
        auto cost = *a - 1.0;
        cost *= cost;

        unsigned i = 1;
        cost += std::inner_product(a+1, b, a, 0.0,
            std::plus<double>(),
            [&i](double x2, double x1) {
                auto t = 2 * x2 * x2 - x1;
                return (++i) * t * t;
            }
        );
        return cost;
    }
    auto domain ( unsigned i ) const -> domain_type { return std::make_pair(-10.0, 10.0); }
    auto extremum ( unsigned i ) const -> double { return 0.0; }
};
