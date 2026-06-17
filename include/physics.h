#include <iostream>
#include <functional>
#include <vector>
#include <chrono>
#include <cstdio>

#include <include/utils.h>

#ifndef PHYSICS_H_INCLUDED
#define PHYSICS_H_INCLUDED


constexpr double PHYS_G = 6.6743015E+04; // km3/(1e24 kg)/s^2, Newtons constant of gravity
constexpr double PHYS_c = 299792.458; // km/s, speed of light
constexpr double PHYS_inv_c2 = 1.0/(PHYS_c*PHYS_c); // 1/c^2


// forward declare Forest_Ruth solver class to be visible here (definition is in solver.h)
template <unsigned N> class Forest_Ruth;


template<unsigned N>
class Hamiltonian {
public:
    Hamiltonian(BodyArray<N> initial_bodies, double t0) :
        bodies_(initial_bodies), t_(t0)
    {
        // initialize p and L
        VectorArray<N> p = p_of_v();
        VectorArray<N> L = L_of_w();
        bodies_.Setp(p);
        bodies_.SetL(L);
    }

    // get the current time
    const double& GetCurrentTime() const {
        return t_;
    }

    // get the current body
    const Body& GetCurrentBody(unsigned idx) const {
        return bodies_[idx];
    }


    VectorArray<N> p_of_v()
    {
        // get the momenta from the velocities
        Compute_lo_metric();

        VectorArray<N> p;
        for (unsigned i = 0; i < N; i++) {
            const double& M = bodies_[i].GetM();
            const Vector& v = bodies_[i].Getv();

            p.row(i) = M*v
                        + PHYS_inv_c2 *(M/2 *v.squaredNorm()*v - 3*M *Phi_[i]*v - 4*M *Theta_[i]);
        }

        return p;
    }


    VectorArray<N> v_of_p()
    {
        // get the velocities from the momenta
        Compute_lo_metric();

        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            const double& M = bodies_[i].GetM();
            const Vector& p = bodies_[i].Getp();

            v.row(i) = p /M
                        + PHYS_inv_c2 *(-1.0/(2*M*M*M) *p.squaredNorm()*p + 3.0/M *Phi_[i]*p + 4.0 *Theta_[i]);
        }

        return v;
    }


    VectorArray<N> L_of_w()
    {
        // get the angular momenta from the angular velocities
        Compute_moment_of_inertia();

        VectorArray<N> L;
        for (unsigned i = 0; i < N; i++) {
            const Vector& w = bodies_[i].Getw();
            L.row(i) = I_[i] *w;
        }

        return L;
    }


    VectorArray<N> w_of_L()
    {
        // get the angular velocities from the angular momenta
        Compute_moment_of_inertia();

        VectorArray<N> w;
        for (unsigned i = 0; i < N; i++) {
            const Vector& L = bodies_[i].GetL();
            w.row(i) = I_inv_[i] *L;
        }

        return w;
    }


    VectorArray<N> dH_kin_dp()
    {
        // derivative of the kinetic Hamiltonian wrt. momentum
        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            const Vector& p = bodies_[i].Getp();
            const double& M = bodies_[i].GetM();
            v.row(i) = p/M
                        + PHYS_inv_c2 *(-1.0/(2*M*M*M) *p.squaredNorm()*p);
        }

        return v;
    }


    VectorArray<N> dH_pot_dx()
    {
        // derivative of the potential Hamiltonian wrt. position
        Compute_lo_metric();
        Compute_grad_Lambda();

        VectorArray<N> NegF = VectorArray<N>::Zero();
        for (unsigned i = 0; i < N; i++) {
            const double& M = bodies_[i].GetM();

            NegF.row(i) = M*grad_Phi_[i]
                            + PHYS_inv_c2 *(M *Phi_[i]*grad_Phi_[i] + M/2 *grad_Lambda_[i]);
        }

        Vector negf;
        Vector r_vec;
        double r;
        Vector e_r;
        double r2_inv;
        double ewxer;
        double ewyer;

        for (unsigned i = 0; i < N; i++) {
            const Vector& x = bodies_[i].Getx();
            const double& Mx = bodies_[i].GetM();
            const double& ax = bodies_[i].Geta();
            const double& J2x = bodies_[i].GetJ2();
            const Vector ewx = bodies_[i].GetAxis();

            for (unsigned j = 0; j < i; j++) {
                const Vector& y = bodies_[j].Getx();
                const double& My = bodies_[j].GetM();
                const double& ay = bodies_[j].Geta();
                const double& J2y = bodies_[j].GetJ2();
                const Vector ewy = bodies_[j].GetAxis();

                r_vec = x - y;
                r = r_vec.norm();
                e_r = r_vec /r;
                r2_inv = 1/(r*r);
                ewxer = ewx.dot(e_r);
                ewyer = ewy.dot(e_r);

                negf = PHYS_G*Mx*My*r2_inv *(
                    -3*J2x*(ax*ax*r2_inv)*((2.5*ewxer*ewxer - 0.5)*e_r - ewxer*ewx)
                    -3*J2y*(ay*ay*r2_inv)*((2.5*ewyer*ewyer - 0.5)*e_r - ewyer*ewy));

                NegF.row(i) += negf;
                NegF.row(j) -= negf;
            }
        }

        return NegF;
    }


    VectorArray<N> dH_kin_dL()
    {
        // derivative of the kinetic Hamiltonian wrt. angular momentum
        return w_of_L();
    }


    VectorArray<N> dH_pot_dorient()
    {
        // derivative of the potential Hamiltonian wrt. orientation
        Compute_moment_of_inertia();

        VectorArray<N> NegTorque = VectorArray<N>::Zero();
        Vector negtorque;
        Vector r_vec;
        double r;
        Vector e_r;
        double r3_inv;
        double ewer;

        for (unsigned i = 0; i < N; i++) {
            const Vector& x = bodies_[i].Getx();
            const double& alpha = bodies_[i].Getorient()(1);
            const double& delta = bodies_[i].Getorient()(2);
            const Vector& L = bodies_[i].GetL();
            const double& Mx = bodies_[i].GetM();
            const double& Ixy = bodies_[i].GetIxy();
            const double& Iz = bodies_[i].GetIz();
            const double& a = bodies_[i].Geta();
            const double& J2 = bodies_[i].GetJ2();
            const Vector ew = bodies_[i].GetAxis();

            const Vector ealpha(-sin(alpha), cos(alpha), 0.0);
            const Vector edelta(-sin(delta)*cos(alpha), -sin(delta)*sin(alpha), cos(delta));

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}
                const Vector& y = bodies_[j].Getx();
                const double& My = bodies_[j].GetM();

                r_vec = y - x;
                r = r_vec.norm();
                e_r = r_vec /r;
                r3_inv = 1/(r*r*r);
                ewer = ew.dot(e_r);

                negtorque(0) = 0.0;
                negtorque(1) = -3*PHYS_G*Mx*My*r3_inv*a*a*J2*ewer*(ealpha.dot(e_r))*cos(delta);
                negtorque(2) = 1.0/2*(I_inv_prime_[i] *L).dot(L)
                            - 3*PHYS_G*Mx*My*r3_inv*a*a*J2*ewer*(edelta.dot(e_r));

                NegTorque.row(i) += negtorque;
            }
        }
        return NegTorque;
    }


    VectorArray<N> dH_pert_dx()
    {
        // derivative of the perturbation Hamiltonian wrt. position
        Compute_lo_metric();

        VectorArray<N> NegF;
        for (unsigned i = 0; i < N; i++) {
            const double& M = bodies_[i].GetM();
            const Vector& p = bodies_[i].Getp();

            NegF.row(i) = PHYS_inv_c2 *(4 *grad_Theta_[i]*p + 3.0/(2*M) *grad_Phi_[i]*p.squaredNorm());
        }

        return NegF;
    }


    VectorArray<N> dH_pert_dp()
    {
        // derivative of the perturbation Hamiltonian wrt. momentum
        Compute_lo_metric();

        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            const double& M = bodies_[i].GetM();
            const Vector& p = bodies_[i].Getp();

            v.row(i) = PHYS_inv_c2 *(4.0 *Theta_[i] + 3.0/M *Phi_[i]*p);
        }

        return v;
    }


private:
    friend class Forest_Ruth<N>;

    BodyArray<N> bodies_;
    double t_;

    heap_array<double, N> Phi_;
    heap_array<Vector, N> Theta_;
    VectorArray<N> pos_Phi_Theta_;

    heap_array<Vector, N> grad_Phi_;
    heap_array<Matrix, N> grad_Theta_;

    heap_array<Vector, N> grad_Lambda_;

    heap_array<Matrix, N> I_;
    heap_array<Matrix, N> I_inv_;
    heap_array<Matrix, N> I_inv_prime_;


    void Compute_lo_metric() {
        // compute the leading order metric perturbations Phi and Theta
        // as well as their gradients

        for (unsigned i = 0; i < N; i++) {
            Phi_[i] = 0.0;
            Theta_[i] = Vector::Zero();
            grad_Phi_[i] = Vector::Zero();
            grad_Theta_[i] = Matrix::Zero();

            const Vector& x0 = bodies_[i].Getx();

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}

                const double& M = bodies_[j].GetM();
                const Vector& x = bodies_[j].Getx();
                const Vector& v = bodies_[j].Getv();

                Vector r_vec = x0 - x;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;

                double Phi_incr = -PHYS_G*M /r;
                Phi_[i] += Phi_incr;
                Theta_[i] += -Phi_incr *v;

                Vector grad_Phi_incr = - Phi_incr /r *e_r;
                grad_Phi_[i] += grad_Phi_incr;
                grad_Theta_[i] += -grad_Phi_incr *v.transpose();
            }
        }
    }


    void Compute_grad_Lambda() {
        // compute the gradient of Lambda

        for (unsigned i = 0; i < N; i++) {
            grad_Lambda_[i] = Vector::Zero();

            const Vector& x0 = bodies_[i].Getx();

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}

                const double& M = bodies_[j].GetM();
                const Vector& x = bodies_[j].Getx();
                const Vector& v = bodies_[j].Getv();

                Vector r_vec = x0 - x;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;

                const Vector& a = -grad_Phi_[j];

                Vector incr_1 = 4*PHYS_G*M *v.squaredNorm() /(r*r) * e_r;
                Vector incr_2 = PHYS_G*M* e_r.dot(v) /(r*r) *(2*v - 3*e_r.dot(v) *e_r);
                Vector incr_3 = -4*PHYS_G*M *Phi_[j] /(r*r) * e_r;
                Vector incr_4 = PHYS_G*M /r *(a - (e_r.dot(a)) *e_r);

                grad_Lambda_[i] += incr_1 + incr_2 + incr_3 + incr_4;
            }
        }
    }


    void Compute_moment_of_inertia()
    {
        // compute the moment of inertia tensor and its inverse
        // as well as the derivative of the inverse wrt. to delta

        for (unsigned i = 0; i < N; i++) {
            const double& delta = bodies_[i].Getorient()(2);
            const double& Ixy = bodies_[i].GetIxy();
            const double& Iz = bodies_[i].GetIz();

            I_[i] << Iz, -Iz*(1.0-sin(delta)), 0.0,
                -Iz*(1.0-sin(delta)), Ixy*cos(delta)*cos(delta) + Iz*(1.0-sin(delta))*(1.0-sin(delta)), 0.0,
                0.0, 0.0,Ixy;

            I_inv_[i] << 1.0/Iz + 1.0/Ixy*cos(delta)*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 1.0/Ixy/(1.0+sin(delta)), 0.0,
                    1.0/Ixy/(1.0+sin(delta)), 1.0/Ixy/cos(delta)/cos(delta), 0.0,
                    0.0, 0.0, 1.0/Ixy;

            I_inv_prime_[i] << -2.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 0.0,
                            -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 2.0/Ixy*sin(delta)/cos(delta)/cos(delta)/cos(delta), 0.0,
                            0.0, 0.0, 0.0;
        }
    }

};









#endif // PHYSICS_H_INCLUDED

