#include <vector>
#include <functional>
#include <cmath>

#include <include/utils.h>
#include <include/physics.h>


#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED


template<unsigned N>
void H_kin_update(double dt, BodyArray<N> &bodies)
{   // do an update of the phase space states according to the kinetic part of the Hamiltonian
    // this must only depend on the momenta
    VectorArray<N> x_dot = dH_kin_dp(bodies);
    VectorArray<N> orient_dot = dH_kin_dL(bodies);

    bodies.Incrementx(x_dot*dt);
    bodies.Incrementorient(orient_dot*dt);
};

template<unsigned N>
void H_pot_update(double dt, BodyArray<N> &bodies)
{   // do an update of the phase space states according to the potential part of the Hamiltonian
    // this must only depend on the coordinated
    VectorArray<N> p_dot = -dH_pot_dx(bodies);
    VectorArray<N> L_dot = -dH_pot_dorient(bodies);

    bodies.Incrementp(p_dot*dt);
    bodies.IncrementL(L_dot*dt);
};

template<unsigned N>
void H_pert_update(double dt, BodyArray<N> &bodies)
{   // do an update of the phase space states according to the perturbation part of the Hamiltonian
    // this can depend on both coordinated and momenta but must be small (i.e. a perturbation)
    VectorArray<N> x0 = bodies.Getx();
    VectorArray<N> p0 = bodies.Getp();

    // do a fixed point iteration for x_mid = (x+x0)/2 starting at x_mid = x0 (and equivalently for p)
    VectorArray<N> x_mid = x0;
    VectorArray<N> p_mid = p0;

    for (unsigned i = 0; i < 2; i++) {
        bodies.Setx(x_mid);
        bodies.Setp(p_mid);

        VectorArray<N> x_dot = dH_pert_dp(bodies);
        VectorArray<N> p_dot = -dH_pert_dx(bodies);

        x_mid = x0 + dt/2 *x_dot;
        p_mid = p0 + dt/2 *p_dot;
    }

    bodies.Setx(2*x_mid - x0);
    bodies.Setp(2*p_mid - p0);
};



// Forest-Ruth 4th order symplectic integrator for N body system
template<size_t N>
class Forest_Ruth {

public:
    // constructor
    Forest_Ruth(
        BodyArray<N> initial_bodies,
        double t0)
        :
        bodies_(initial_bodies),
        t_(t0)
    {
        // initialize p and L
        VectorArray<N> p = p_of_v(bodies_);
        VectorArray<N> L = L_of_w(bodies_);
        bodies_.Setp(p);
        bodies_.SetL(L);
    }

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

        // keep v and w up to date
        VectorArray<N> v = v_of_p(bodies_);
        VectorArray<N> w = w_of_L(bodies_);
        bodies_.Setv(v);
        bodies_.Setw(w);
    }


private:
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
        H_pot_update(dt/2, bodies);

        // step 2: update using H_kin by half a step
        H_kin_update(dt/2, bodies);

        // step 3: update using H_pert by a full step
        H_pert_update(dt, bodies);

        // step 4: update unsing H_kin by another half a step
        H_kin_update(dt/2, bodies);

        // step 5: update using H_pot by another half a step
        H_pot_update(dt/2, bodies);
    }

    // compute a step of the Forest-Ruth algorithm
    void ComputeForestRuthStep(double dt, BodyArray<N> &bodies) {

        ComputeLeapFrogStep(c1*dt, bodies);
        ComputeLeapFrogStep(c2*dt, bodies);
        ComputeLeapFrogStep(c1*dt, bodies);
    }

};



#endif //SOLVER_H_INCLUDED

