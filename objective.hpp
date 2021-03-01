#ifndef HPP_OBJECTIVE
#define HPP_OBJECTIVE

#include <algorithm>
#include <stdexcept>
#include <utility>

//    using typename objective<FwdIter>::domain_type;
//    using param_type = std::vector<double>::const_iterator;

class objective
{
public:
    using domain_type = std::pair<double,double>;
    using param_type = std::vector<double>::const_iterator;
    using result_type = double;

    virtual auto
    operator() ( param_type a, param_type b ) const -> result_type = 0;

    virtual auto domain() const -> domain_type = 0;
};

class sphere : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
    {
        auto cost = std::inner_product(a, b, a, 0.0);
        return cost;
    }

    auto domain() const -> domain_type { return std::make_pair(-5.12, 5.12); }
};

class rosenbrock : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
    {
        std::vector<double> tmp;
        std::adjacent_difference(a, b, std::back_inserter(tmp),
            [](double x2, double x1) {
                double t1 = x2 - x1 * x1;
                double t2 = 1 - x1;
                return 100.0 * t1 * t1 + t2 * t2;
            }
        );
        double cost = std::accumulate(tmp.begin()+1, tmp.end(), 0.0);
        return cost;
    }

    auto domain() const -> domain_type { return std::make_pair(-5.0, 10.0); }
};

class rastrigin : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
    {
        double cost = std::accumulate(a, b, 0.0,
            [](double sum, double x) {
                static double const TWOPI = 8 * std::atan(1.0);
                return sum + x * x - 10 * std::cos(TWOPI * x);
            }
        );            
        unsigned n = distance(a,b);
        return 10 * n + cost;
    }
    auto domain() const -> domain_type { return std::make_pair(-5.12, 5.12); }
};

class griewangk : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
    {
        double cost1 = std::accumulate(a, b, 0.0,
            [](double sum, double x) { return sum + x * x  / 4000.0; }
        );
        unsigned i = 0;
        double cost2 = std::accumulate(a, b, 1.0,
            [&i](double prod, double x) {
                return prod * std::cos(x / std::sqrt(++i));
            }
        );
        return cost1 - cost2 + 1;
    }
    auto domain() const -> domain_type { return std::make_pair(-600.0, 600.0); }
};

class shaffer_f6 : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
    {
        auto x = *(a++), y = *a;
        auto h = x * x + y * y;
        auto denom = 1 + 0.001 * h;
        auto numer = std::sin(std::sqrt(h));  
        double cost = 0.5 + (numer * numer - 0.5) / (denom * denom);
        return cost;
    }
    //auto domain() const -> domain_type { return std::make_pair(-600.0, 600.0); }
};

class ackley : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
    {
        auto n = std::distance(a, b);
        constexpr auto TWOPI = 8.0 * std::atan(1.0);
        double s1 = std::accumulate(a, b, 0.0,
            [](double sum, double x) { return sum + x * x; }
        );
        double s2 = std::accumulate(a, b, 0.0,
            [](double sum, double x) { return sum + std::cos(TWOPI * x); }
        );
        double cost = 20.0 + std::exp(1.0)
                       - 20.0 * std::exp(-0.2 * std::sqrt(1.0 / n * s1))
                       - std::exp(1.0 / n * s2);
        return cost;
    }
    auto domain() const -> domain_type { return std::make_pair(-15.0, 30.0); }
};

class beale : public objective
{
public:
    auto operator() ( param_type a, param_type b ) const -> result_type
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
    auto domain() const -> domain_type { return std::make_pair(-4.5, 4.5); }
};

#endif
