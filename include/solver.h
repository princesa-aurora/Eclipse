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
        std::function<void(double, BodyArray<N>&)> H_kin_update,
        std::function<void(double, BodyArray<N>&)> H_pot_update,
        std::function<void(double, BodyArray<N>&)> H_pert_update,
        BodyArray<N> initial,
        double t0)
        :
        H_kin_update_(H_kin_update),
        H_pot_update_(H_pot_update),
        H_pert_update_(H_pert_update),
        bodies_(initial),
        t_(t0)
    {}

    // destructor
    ~Forest_Ruth() = default;

    // get the last time
    const double& GetCurrentTime() const {
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
    // function that updates the phase space states according to the kinetic Hamiltonian
    const std::function<void(double, BodyArray<N>&)> H_kin_update_;
    // function that updates the phase space states according to the potential Hamiltonian
    const std::function<void(double, BodyArray<N>&)> H_pot_update_;
    // function that updates the phase space states according to the perturbation Hamiltonian
    const std::function<void(double, BodyArray<N>&)> H_pert_update_;

    // current time
    double t_;
    // current bodies
    BodyArray<N> bodies_;

    // algorithm constants
    static constexpr double c1 = 1/(2-cbrt(2));
    static constexpr double c2 = 1-2*c1;

    // compute a step of the Leapfrog algorithm
    void ComputeLeapFrogStep(double dt, BodyArray<N> &bodies) {

        // step 1: update using H_pot by half a step
        H_pot_update_(dt/2, bodies);

        // step 2: update using H_kin by half a step
        H_kin_update_(dt/2, bodies);

        // step 3: update using H_pert by a full step
        H_pert_update_(dt, bodies);

        // step 4: update unsing H_kin by another half a step
        H_kin_update_(dt/2, bodies);

        // step 5: update using H_pot by another half a step
        H_pot_update_(dt/2, bodies);
    }

    // compute a step of the Forest-Ruth algorithm
    void ComputeForestRuthStep(double dt, BodyArray<N> &bodies) {

        ComputeLeapFrogStep(c1*dt, bodies);
        ComputeLeapFrogStep(c2*dt, bodies);
        ComputeLeapFrogStep(c1*dt, bodies);
    }

};



#endif //SOLVER_H_INCLUDED

