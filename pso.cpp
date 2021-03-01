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

class swarm;

class particle : protected std::vector<double>
{
protected:
    using super = std::vector<double>;
    using param = super::const_iterator;
    using solution = std::pair<super::value_type,super>;
    super velocity;
public:
    //using super::vector; // this should replace the constructors below
    particle() : super() {}
    particle(unsigned n)
        : super(n), velocity(n,0),
          local_best(std::numeric_limits<double>::infinity(),*this)
        {}
    using super::operator[];
    using super::begin;
    using super::cbegin;
    using super::end;
    using super::cend;
    using super::size;
    using super::const_iterator;
    friend bool operator< ( solution const& a, solution const& b );
    friend bool operator< ( double a, solution const& b );
    friend bool operator< ( solution const& a, double b );
private:
    friend swarm;
    solution local_best;
};

bool operator< ( particle::solution const& a, particle::solution const& b ) { return a.first < b.first; }
bool operator< ( double a, particle::solution const& b ) { return a < b.first; }
bool operator< ( particle::solution const& a, double b ) { return a.first < b; }

struct objective
{
    using param = particle::const_iterator;
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
        auto cost = std::inner_product(a+1, b, a, 0.0,
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

class swarm : std::vector<particle>
{
    using super = std::vector<particle>;
public:
    struct param_type
    {
        /* Population size */
        double n;
        /* Cognitive trust parameter */
        double c1;
        /* Social trust parameter */
        double c2;
        /* Current inertia */
        double w;
        /* Initial inertia */
        double w0;
        /* Inertial decay */
        double wd;
        /* Velocity fraction */
        double k;
        /* Velocity decay */
        double vd;
        /* Decay delay in iterations */
        long d;
    };

    explicit swarm ( int d, objective* f = new sphere(),
                      param_type p = {20,2.0,2.0,1.0,1.0,0.95,0.5,0.95,200} )
        : super(p.n, d), param(p), f(f), leader(begin()),
          vmax(d, 0.0), rng(std::random_device()())
    {
        auto i = 0L;
        std::generate(vmax.begin(), vmax.end(), [&i,this](){
            auto bounds = this->f->domain(i++);
            auto range = bounds.second - bounds.first;
            return range * param.k;
        });
    }

    particle::solution best_solution() { return leader->local_best; }

    void operator() ()
    {
        std::cerr << param.n << ' '
                  << param.c1 << ' '
                  << param.c2 << ' '
                  << param.w << ' '
                  << param.w0 << ' '
                  << param.wd << ' '
                  << param.k << ' '
                  << param.vd << ' '
                  << param.d << std::endl;
        initialize();

        auto k = 0LL, t = 0LL;
        do {
            for ( auto i = begin(); i != end(); ++i, ++k ) {
                if ( update(i) ) {
                    t = 0;
                } else {
                    ++t;
                }
                if (t == param.d) {
                    t = 0;
                    param.w  *= param.wd;
                    for ( auto& vm : vmax ) {
                        vm *= param.vd;
                    }
                }

                if (leader->local_best.first < 0.1 || k > 640000) {
                    std::cerr << k << std::endl;
                    return;
                }
            }
        } while (true);
    }

private:
    bool update ( iterator i )
    {
        // compute velocity
        {
            auto xi = i->cbegin();
            auto vm = vmax.cbegin();
            auto p = i->local_best.second.cbegin();
            auto g = leader->local_best.second.cbegin();
            std::transform(i->velocity.cbegin(), i->velocity.cend(),
                           i->velocity.begin(),
                           [&,this](double vprev){
                auto r1 = std::generate_canonical<double,16>(rng);
                auto r2 = std::generate_canonical<double,16>(rng);
                auto x = *xi++;
                auto vmax = *vm++;
                auto v = vprev * param.w
                    + r1 * param.c1 * (*p++ - x)
                    + r2 * param.c2 * (*g++ - x);
                return std::max( std::min( v, vmax ), -vmax );
            });
        }

        // update position
        std::transform(i->cbegin(), i->cend(), i->velocity.cbegin(),
                       i->begin(), std::plus<double>());
        // compute cost
        auto cost = (*f)(i->begin(), i->end());
        // update personal best
        if ( cost < i->local_best ) {
            i->local_best.second.assign(i->cbegin(),i->cend());
            i->local_best.first = cost;
            // update global best
            if ( cost < leader->local_best ) {
                leader = i;
            }
            if (leader == i) { return true; }
        }
        return false;
    }

    void initialize()
    {
        randomize();
        for ( auto i = begin(); i != end(); ++i) {
            // compute cost
            auto cost = (*f)(i->begin(), i->end());
            // update personal best
            i->local_best.second.assign(i->cbegin(),i->cend());
            i->local_best.first = cost;
            if (cost < leader->local_best) {
                leader = i;
            }
        }
    }

    void randomize()
    {
        std::uniform_real_distribution<> dis;
        for ( auto& p : *this )
        {
            auto i = 0U;
            std::generate(p.begin(), p.end(), [&i,&dis,this](){
                dis.param(std::uniform_real_distribution<>::param_type(
                    f->domain(i).first, f->domain(i).second
                ));
                ++i;
                return dis(rng);
            });
        }
        for ( auto& p : *this )
        {
            auto vm = vmax.cbegin();
            std::generate(p.velocity.begin(), p.velocity.end(), [&,this](){
                dis.param(std::uniform_real_distribution<>::param_type( -*vm, *vm ));
                ++vm;
                return dis(rng);
            });
        }
    }

    param_type param;
    std::unique_ptr<objective> f;
    super::iterator leader;
    std::vector<double> vmax;
    std::mt19937 rng;
};

int main (int argc, char** argv)
{
    /*
    int const N = atoi(argv[1]);
    double w0, c1, c2;
    std::istringstream arg(argv[2]);
    arg >> w0;
    arg.str(argv[3]);
    arg >> c1;
    arg.str(argv[4]);
    arg >> c2;
    */

    swarm s { 64, new griewangk() };
    s();

    //double bestcost = s.best_solution().first;
    //auto const& best = s.best_solution().second;
    //std::cout << bestcost << std::endl;
    //std::copy(best.cbegin(), best.cend(), std::ostream_iterator<double>(std::cout," "));
    //std::cout << std::endl;

    return 0;
}
