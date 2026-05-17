#include <vector>
#include <functional>
#include <cmath>

#include <include/utils.h>


#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED


// Forest-Ruth 4th order symplectic integrator for N body system
template<size_t N>
class Forest_Ruth {

public:
    // constructor
    Forest_Ruth(
        std::function<VectorArray<N>(const BodyArray<N>&)> x_dot,
        std::function<VectorArray<N>(const BodyArray<N>&)> p_dot,
        std::function<VectorArray<N>(const BodyArray<N>&)> orient_dot,
        std::function<VectorArray<N>(const BodyArray<N>&)> L_dot,
        BodyArray<N> initial,
        double t0)
        :
        x_dot_(x_dot),
        p_dot_(p_dot),
        orient_dot_(orient_dot),
        L_dot_(L_dot),
        bodies_(initial),
        t_(t0)
    {}

    // destructor
    ~Forest_Ruth() = default;

    // get the last time
    double GetCurrentTime() const {
        return t_;
    }

    // get the current bodies
    const Body& GetCurrentBody(unsigned idx) const {
        return bodies_[idx];
    }

    // make a step of the Forest Ruth algorithm
    void MakeStep(double dt) {

        t_ += dt;
        ComputeForestRuthStep(dt, bodies_);
    }


private:
    // function that gives x_dot(x, p)
    std::function<VectorArray<N>(const BodyArray<N>&)> x_dot_;
    // function that gives p_dot(x, p)
    std::function<VectorArray<N>(const BodyArray<N>&)> p_dot_;
    // function that gives orient_dot(x, p)
    std::function<VectorArray<N>(const BodyArray<N>&)> orient_dot_;
    // function that gives L_dot(x, p)
    std::function<VectorArray<N>(const BodyArray<N>&)> L_dot_;

    // vector that stores the times for which values were calculated
    std::vector<double> T_;
    // vector that stores the computed bodies
    std::vector<BodyArray<N>> Bodies_;

    // current time
    double t_;
    // current bodies
    BodyArray<N> bodies_;

    // algorithm constants
    static constexpr double c1 = 1/(2-cbrt(2));
    static constexpr double c2 = 1-2*c1;

    // compute a step of the Leapfrog algorithm
    void ComputeLeapFrogStep(double dt, BodyArray<N> &bodies) {

        // step 1: update momenta by half a step
        bodies.Incrementp(dt/2 * p_dot_(bodies));
        bodies.IncrementL(dt/2 * L_dot_(bodies));

        // step 2: update positions by a full step
        bodies.Incrementx(dt * x_dot_(bodies));
        bodies.Incrementorient(dt * orient_dot_(bodies));

        // step 3: update momenta by another half a step
        bodies.Incrementp(dt/2 * p_dot_(bodies));
        bodies.IncrementL(dt/2 * L_dot_(bodies));
    }

    // compute a step of the Forest-Ruth algorithm
    void ComputeForestRuthStep(double dt, BodyArray<N> &bodies) {

        ComputeLeapFrogStep(c1*dt, bodies);
        ComputeLeapFrogStep(c2*dt, bodies);
        ComputeLeapFrogStep(c1*dt, bodies);
    }

};



#endif //SOLVER_H_INCLUDED

