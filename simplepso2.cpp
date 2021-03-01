#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <random>
#include <typeinfo>
#include <vector>

long P = 20;
long N = 64;
long kappa;
long kmax;
double c1;
double c2;
double w0;
double vmax;
double w;
double vd;
double wd;
long k = 0;
long t = 0;
long d = 200;
std::mt19937 rng;
std::uniform_real_distribution<> dis;
std::vector<std::vector<double>> x;
std::vector<std::vector<double>> v;
std::vector<std::vector<double>> p;
std::vector<double> f;
std::vector<double> g;
double fg;

double cost( std::vector<double> const& );
void initialize();
void optimize();
void report();

using position = std::vector<double>;
using velocity = std::vector<double>;

position& operator+= (position& l, velocity const& r);

int main()
{
    initialize();
    optimize();
    report();

    return 0;
}

void initialize()
{
    // a. Set constants κ , c1, c2, kmax,, w0, νd, wd, and d
    k = 0;
    c1 = 2.0;
    c2 = 2.0;
    kmax = 640000;
    vmax = (600.0 - -600.0) * 0.5;
    w = 1.0;
    vd = 1.0 - 0.05;
    wd = 1.0 - 0.05;
    d = 200;

    // b. Set counters k=0, t=0. Set random number seed
    k = 0;
    t = 0;
    rng.seed(std::random_device()());

    // c. Randomly initialize particle positions  in  for i=1,... , p
    x.resize(P);
    dis.param(std::uniform_real_distribution<>::param_type(-600.0, 600.0));
    for ( auto& xk : x )
    {
        std::generate_n(std::back_inserter(xk), N, [](){
            return dis(rng);
        });
    }

    // d. Randomly initialize particle velocities  for i=1,..., p
    v.resize(P);
    dis.param(std::uniform_real_distribution<>::param_type(-vmax, vmax));
    for ( auto& vk : v )
    {
        std::generate_n(std::back_inserter(vk), N, [](){
            return dis(rng);
        });
    }
    
    // e. Evaluate cost function values using design space coordinates for i=1,..., p
    std::transform(x.cbegin(), x.cend(), std::back_inserter(f), cost);

    // f. Set and for i=1,..., p
    std::copy(x.cbegin(), x.cend(), std::back_inserter(p));

    // g. Set to best and g0 to corresponding
    auto best = std::min_element(f.cbegin(), f.cend());
    fg = *best;
    auto i = x.cbegin();
    std::advance(i, std::distance(f.cbegin(), best));
    g = *i;

    std::cout << fg << ": ";
    std::copy(f.cbegin(), f.cend(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
}

double cost( std::vector<double> const& xk )
{
    auto a = xk.cbegin();
    auto b = xk.cend();
    auto cost1 = std::accumulate(a, b, 0.0,
        [](double sum, double x) { return sum + x * x  / 4000.0; }
    );
    auto i = 0.0;
    auto cost2 = std::accumulate(a, b, 1.0,
        [&i](double prod, double x) {
            return prod * std::cos(x / std::sqrt(++i));
        }
    );
    ++k;
    return cost1 - cost2 + 1.0;
}

void optimize()
{
    dis.param(std::uniform_real_distribution<>::param_type(0.0, 1.0));

    while (true)
    {
        for ( auto i = 0U; i < x.size(); ++i )
        {
            auto& xi = x[i];
            auto& vi = v[i];
            auto& pi = p[i];

            // a. Update particle velocity vectors  using Eq. (2)
            auto pk = pi.cbegin(), xj = xi.cbegin(), gk = g.cbegin();
            //auto r1 = dis(rng), r2 = dis(rng);
            std::transform(vi.cbegin(), vi.cend(), vi.begin(),
                [&pk,&xj,&gk](double vk){
                    // b. If  for any component, then set that component to its maximum allowable value
                    auto r1 = dis(rng), r2 = dis(rng);
                    auto xk = *xj++;
                    return std::max(-vmax, std::min(vmax, w * vk + c1 * r1 * (*pk++ - xk) + c2 * r2 * (*gk++ - xk)));
                }
            );
    
            // c. Update particle position vectors  using Eq. (1)
            std::transform(xi.cbegin(), xi.cend(), vi.cbegin(), xi.begin(), std::plus<double>());
        }
    
        // d. Evaluate cost function values  using design space coordinates  for i=1,..., p
        // e. If , then ,  for i=1,..., p
        auto pk = p.begin();
        std::transform(x.cbegin(), x.cend(), f.cbegin(), f.begin(),
            [&pk](std::vector<double> const& xk, double fj){
                auto fk = cost(xk);
                if (fk < fj) {
                    std::copy(xk.cbegin(), xk.cend(), pk->begin());
                    fj = fk;
                }
                ++pk;
                return fj;
        });

        // f. If  then , for i=1,..., p
        std::cout << fg << ": ";
        std::copy(f.cbegin(), f.cend(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;

        auto best = std::min_element(f.cbegin(), f.cend());
        auto i = x.cbegin();
        std::advance(i, std::distance(f.cbegin(), best));
        if ( *best < fg ) {
            fg = *best;
            g = *i;
            // g. If  was improved in (e), then reset t=0, otherwise increment t
            t = 0;
            //std::cerr << '!';
        } else {
            ++t;
            //std::cerr << '.';
        }

        // h. If the maximum number of function evaluations is exceeded, then go to 3
        if ( k > kmax || fg < 0.1 ) break;

        // i. If t=d, then multiply wk+1 by (1-wd) and  by (1 -νd)
        if (t == d) {
            //std::cerr << '?';
            w *= wd;
            vmax *= vd;
        }
    }
}

void report()
{
    std::cout << "k: " << k << std::endl;
    std::cout << "fg: " << fg << std::endl;
    std::cout << "g: ";
    std::copy( g.begin(), g.end(), std::ostream_iterator<double>(std::cout," ") );
    std::cout << std::endl;
}

velocity& operator*= (velocity& l, double r)
{
    using namespace std::placeholders;
    std::transform ( l.cbegin(), l.cend(), l.begin(),
                     std::bind(std::multiplies<double>(), _1, r) );
    return l;
}

position& operator+= (position& l, velocity const& r)
{
    std::transform ( l.cbegin(), l.cend(), r.cbegin(), l.begin(), std::plus<double>() );
    return l;
}
