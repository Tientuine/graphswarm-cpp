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
    particle(unsigned n) : super(n), velocity(n,0) {}
    //particle(std::initializer_list<double> x) : super(x) {}
    using super::operator[];
    using super::begin;
    using super::cbegin;
    using super::end;
    using super::cend;
    using super::size;
    using super::const_iterator;
    friend bool operator< ( solution a, solution b );
    friend bool operator< ( double a, solution b );
    friend bool operator< ( solution a, double b );
};

bool operator< ( particle::solution a, particle::solution b ) { return a.first < b.first; }
bool operator< ( double a, particle::solution b ) { return a < b.first; }
bool operator< ( particle::solution a, double b ) { return a.first < b; }

struct objective
{
    using param = particle::const_iterator;
    using domain_type = std::pair<double,double>;

    virtual auto operator() ( param a, param b ) const -> double = 0;
    virtual auto domain ( unsigned i ) const -> domain_type = 0;
    virtual auto extremum ( unsigned i ) const -> double = 0;
    //virtual auto tolerance () const -> double = 0;
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
    auto tolerance () const -> double { return 0.01; }
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

class swarmer
    : public runnable, particle, public std::enable_shared_from_this<swarmer> 
{
public:
    using swarm = std::forward_list<std::shared_ptr<swarmer>>;

    static long long update_count;

    explicit swarmer ( int n = 1, objective* f = new sphere() )
        : particle(n),
          leader(std::shared_ptr<swarmer>(this)),
          neighbors(new swarm(1,shared_from_this())),
          cost_function(f),
          local_best(std::numeric_limits<double>::infinity(),*this),
          rng(std::random_device()())
        {
            auto i = 0U;
            std::uniform_real_distribution<> dis;
            std::generate(begin(), end(), [&i,&dis,this](){
                dis.param(std::uniform_real_distribution<>::param_type(
                    cost_function->domain(i).first, cost_function->domain(i).second
                ));
                ++i;
                return dis(this->rng);
            });
            std::generate(velocity.begin(), velocity.end(), [&i,&dis,this](){
                auto range =
                    std::abs(cost_function->domain(i).second - cost_function->domain(i).first);
                dis.param(std::uniform_real_distribution<>::param_type( -range, range ));
                ++i;
                return dis(this->rng);
            });

        }

    swarmer ( swarmer & s )
        : particle(s.size()),
          leader(s.shared_from_this()),
          cost_function(s.cost_function),
          local_best(std::numeric_limits<double>::infinity(),*this),
          rng(std::random_device()())
        {
            auto i = 0U;
            std::uniform_real_distribution<> dis;
            std::generate(begin(), end(), [&i,&dis,this](){
                dis.param(std::uniform_real_distribution<>::param_type(cost_function->domain(i).first, cost_function->domain(i).second));
                ++i;
                return dis(this->rng);
            });
            std::generate(velocity.begin(), velocity.end(), [&i,&dis,this](){
                auto range =
                    std::abs(cost_function->domain(i).second - cost_function->domain(i).first);
                dis.param(std::uniform_real_distribution<>::param_type( -range, range ));
                ++i;
                return dis(this->rng);
            });
            leader->neighbors->push_front(std::shared_ptr<swarmer>(this));
        }

    void start()
    {
        std::lock_guard<std::mutex> lock(leader_mutex);
        for ( auto s : *(leader->neighbors) ) { s->run(); }
    }

    void watch() { for ( auto s : *(leader->neighbors) ) { s->join(); } }

    solution best_solution() { return leader->local_best; }

private:
    static bool const DEBUG = true;

    void operator() ()
    {
        while (leader->local_best.first > 0.1)
            { ++update_count; update(); }
    }

    void update()
    {
        // compute velocity
        auto pcurr = cbegin();
        {
            std::lock_guard<std::mutex> lock(leader_mutex);
            auto pbest = local_best.second.cbegin();
            auto gbest = leader->local_best.second.cbegin();
            std::transform(velocity.cbegin(), velocity.cend(),
                           velocity.begin(),
                           [&,this](double v){
                auto prand = std::generate_canonical<double,16>(this->rng);
                auto grand = std::generate_canonical<double,16>(this->rng);

                auto newv = v * this->INERTIA
                        + prand * this->P_AFFINITY * (*(pbest++) - *(pcurr))
                        + grand * this->G_AFFINITY * (*(gbest++) - *(pcurr++));
                return newv;
            });
        }
        std::this_thread::yield();

        // update position
        std::transform(cbegin(), cend(), velocity.cbegin(),
                       begin(), std::plus<double>());
        // compute cost
        auto cost = (*cost_function)(begin(), end());
        // update personal best
        if ( cost < local_best ) {
            local_best.second.assign(cbegin(),cend());
            local_best.first = cost;
        }
        // update global best
        if ( local_best < leader->local_best ) {
            lead();
            std::this_thread::yield();
         }
    }

    void lead()
    { 
        std::lock_guard<std::mutex> lock(leader_mutex);

        neighbors = std::move(leader->neighbors);

        for ( auto s : *neighbors )
            { if (s.get() != this) { s->follow(*this); } }
        follow(*this);

        //std::cerr << leader->local_best.first << std::endl;
    }

    void follow( swarmer & s ) { leader = s.shared_from_this(); }

    std::shared_ptr<swarmer> leader;
    std::unique_ptr<swarm> neighbors;
    std::shared_ptr<objective> cost_function;
    solution local_best;
    std::mt19937 rng;

    static std::mutex leader_mutex;
public:
    static double INERTIA;
    static double P_AFFINITY;
    static double G_AFFINITY;
};

double swarmer::INERTIA = 0.8;
double swarmer::P_AFFINITY = 2;
double swarmer::G_AFFINITY = 2;

std::mutex swarmer::leader_mutex;

long long swarmer::update_count {0};

int main (int argc, char** argv)
{
    int const N = atoi(argv[1]);

    std::istringstream arg(argv[2]);
    arg >> swarmer::INERTIA;
    arg.str(argv[3]);
    arg >> swarmer::P_AFFINITY;
    arg.str(argv[4]);
    arg >> swarmer::G_AFFINITY;

    //std::shared_ptr<swarmer> s {new swarmer(2, new shaffer_f6())};
    std::shared_ptr<swarmer> s {new swarmer(64, new griewangk())};
    //std::shared_ptr<swarmer> s {new swarmer(10, new rosenbrock())};
    for ( auto i = 1; i < N; ++i ) {
        new swarmer(*s);
    }
    try {
        s->start();
        s->watch();
    } catch (std::system_error e) {
        std::cerr << e.code() << std::endl;
    } catch (std::exception e) {
        std::cerr << e.what() << std::endl;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
    }

    double bestcost = s->best_solution().first;
    auto const& best = s->best_solution().second;
    std::cout << swarmer::update_count << std::endl;
    std::cout << bestcost << std::endl;
    std::copy(best.cbegin(), best.cend(), std::ostream_iterator<double>(std::cout," "));
    std::cout << std::endl;
    std::cout << swarmer::INERTIA << std::endl;

    return 0;
}
