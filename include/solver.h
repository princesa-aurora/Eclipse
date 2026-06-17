#include <vector>
#include <functional>
#include <cmath>

#include <include/utils.h>
#include <include/physics.h>


#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED



// Forest-Ruth 4th order symplectic integrator for N body system
template<size_t N>
class Forest_Ruth {

public:
    // constructor
    Forest_Ruth(
        Hamiltonian<N> H)
        :
        H_(H)
    {
        // initialize p and L
        VectorArray<N> p = H_.p_of_v();
        VectorArray<N> L = H_.L_of_w();
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
        ComputeForestRuthStep(dt);

        // keep v and w up to date
        VectorArray<N> v = H_.v_of_p();
        VectorArray<N> w = H_.w_of_L();
        bodies_.Setv(v);
        bodies_.Setw(w);
    }


private:
    // Hamiltonian instance
    Hamiltonian<N> H_;
    // current time
    double& t_ = H_.Gett();
    // current bodies
    BodyArray<N>& bodies_ = H_.Getbodies();

    // algorithm constants
    static constexpr double c1 = 1/(2-cbrt(2));
    static constexpr double c2 = 1-2*c1;

    // compute a step of the Leapfrog algorithm
    void ComputeLeapFrogStep(double dt) {

        // step 1: update using H_pot by half a step
        H_pot_update(dt/2);

        // step 2: update using H_kin by half a step
        H_kin_update(dt/2);

        // step 3: update using H_pert by a full step
        H_pert_update(dt);

        // step 4: update unsing H_kin by another half a step
        H_kin_update(dt/2);

        // step 5: update using H_pot by another half a step
        H_pot_update(dt/2);
    }

    // compute a step of the Forest-Ruth algorithm
    void ComputeForestRuthStep(double dt) {

        ComputeLeapFrogStep(c1*dt);
        ComputeLeapFrogStep(c2*dt);
        ComputeLeapFrogStep(c1*dt);
    }





    void H_kin_update(double dt)
    {   // do an update of the phase space states according to the kinetic part of the Hamiltonian
        // this must only depend on the momenta
        VectorArray<N> x_dot = H_.dH_kin_dp();
        VectorArray<N> orient_dot = H_.dH_kin_dL();

        bodies_.Incrementx(x_dot*dt);
        bodies_.Incrementorient(orient_dot*dt);
    }


    void H_pot_update(double dt)
    {   // do an update of the phase space states according to the potential part of the Hamiltonian
        // this must only depend on the coordinated
        VectorArray<N> p_dot = -H_.dH_pot_dx();
        VectorArray<N> L_dot = -H_.dH_pot_dorient();

        bodies_.Incrementp(p_dot*dt);
        bodies_.IncrementL(L_dot*dt);
    }


    void H_pert_update(double dt)
    {   // do an update of the phase space states according to the perturbation part of the Hamiltonian
        // this can depend on both coordinated and momenta but must be small (i.e. a perturbation)
        VectorArray<N> x0 = bodies_.Getx();
        VectorArray<N> p0 = bodies_.Getp();

        // do a fixed point iteration for x_mid = (x+x0)/2 starting at x_mid = x0 (and equivalently for p)
        VectorArray<N> x_mid = x0;
        VectorArray<N> p_mid = p0;

        for (unsigned i = 0; i < 2; i++) {
            bodies_.Setx(x_mid);
            bodies_.Setp(p_mid);

            VectorArray<N> x_dot = H_.dH_pert_dp();
            VectorArray<N> p_dot = -H_.dH_pert_dx();

            x_mid = x0 + dt/2 *x_dot;
            p_mid = p0 + dt/2 *p_dot;
        }

        bodies_.Setx(2*x_mid - x0);
        bodies_.Setp(2*p_mid - p0);
    }




};



#endif //SOLVER_H_INCLUDED

