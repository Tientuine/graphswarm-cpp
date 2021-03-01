#ifndef HPP_PARTICLE
#define HPP_PARTICLE

#include <algorithm>
#include <forward_list>
#include <functional>
#include <iostream>
#include <memory> 
#include <mutex>
#include <numeric>
#include <random>
#include <vector>
#include "runnables.hpp"

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
    particle(std::initializer_list<double> x) : super(x) {}
    using super::operator[];
    using super::begin;
    using super::cbegin;
    using super::end;
    using super::cend;
    using super::size;
    friend bool operator< ( solution a, solution b );
};

bool operator< ( particle::solution a, particle::solution b )
{ return a.first < b.first; }


class swarmer
    : public std::enable_shared_from_this<swarmer>,
      public runnable,
      particle,
      std::recursive_mutex
{
    using swarm = std::forward_list<std::shared_ptr<swarmer>>;
    using cost_function = std::function<double(param,param)>;
public:
    explicit swarmer ( cost_function f =
        std::bind(std::accumulate<std::vector<double>::const_iterator,
                                double,
                                std::plus<double>
                                >,
                                std::placeholders::_1,
                                std::placeholders::_2,
                                0,std::plus<double>())
     )
        : objective(f),
          leader(this),
          neighbors(new swarm(1,std::shared_ptr<swarmer>(this))),
          local_best(std::numeric_limits<double>::infinity(),*this)
    {}

    explicit swarmer ( swarmer const& p )
        : objective(p.objective),
          leader(const_cast<swarmer&>(p).shared_from_this()),
          local_best(std::numeric_limits<double>::infinity(),*this)
    { leader->neighbors->push_front(std::shared_ptr<swarmer>(this)); }

    void start_swarming() { for ( auto s : *(leader->neighbors) ) { s->run(); } }
    void watch() { for ( auto s : *(leader->neighbors) ) { s->join(); } }
    using runnable::join;

protected:
    void operator() () { while (true) { update(); } }

    void update()
    {
        std::cerr << get_id() << std::endl;
        auto cost = objective(begin(), end());
        if ( cost < local_best.first )
            { local_best = solution(cost,*this); }
        if ( local_best < leader->local_best )
            { promote(); }
    }

    void promote ()
    {
        leader->lock();
        swarmer& former_leader = *leader;
        for ( auto s : *(former_leader.neighbors) )
            { s->leader = shared_from_this(); }
        neighbors = std::move(former_leader.neighbors);
        former_leader.unlock();
    }

    void demote()
    {
        /*
        auto pbest = local_best.second.cbegin();
        auto gbest = global_best().second.cbegin();
        std::transform(velocity.cbegin(), velocity.cend(),
                       velocity.begin(),
                       [&,this](double v){
            auto prand = std::generate_canonical<double,16>(this->rng);
            auto grand = std::generate_canonical<double,16>(this->rng);
            return v * this->INERTIA
                    + prand * this->P_AFFINITY * (*(pbest.operator++(0)))
                    + grand * this->G_AFFINITY * (*(gbest++));
        });
        std::transform(cbegin(), cend(), velocity.cbegin(),
                       begin(), std::plus<double>());
        */
    }

    cost_function objective;
    std::shared_ptr<swarmer> leader;
    std::unique_ptr<swarm> neighbors;
    solution local_best;
    std::mt19937 rng;

    constexpr static double const INERTIA = 0.86;
    constexpr static double const P_AFFINITY = 0.45;
    constexpr static double const G_AFFINITY = 0.25;
};

#endif //HPP_PARTICLE
