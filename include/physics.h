#include <iostream>
#include <functional>
#include <vector>
#include <chrono>
#include <cstdio>

#include <include/utils.h>

#ifndef PHYSICS_H_INCLUDED
#define PHYSICS_H_INCLUDED



// forward declare Forest_Ruth solver class to be visible here (definition is in solver.h)
template <unsigned N> class Forest_Ruth;


template<unsigned N>
class Hamiltonian {
public:
    Hamiltonian(BodyArray<N> initial_bodies, double t0) :
        bodies_(initial_bodies), t_(t0)
    {
        VectorArray<N> p_0th;
        QuaternionArray<N> L_0th;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& M = body.GetM();
            const Vector& v = body.Getv();
            const double& Ixy = body.GetIxy();
            const Quaternion& w = body.Getw();
            p_0th.row(i) = M*v;
            L_0th.row(i) = 4*Ixy*w;

            lambda0_[i] = 1.0/(8*Ixy)* L_0th.row(i).squaredNorm();
        }
        // seed p and L with the 0th order approximations to use for perturbation contributions that depend on v(p) and w(L)
        bodies_.Setp(p_0th);
        bodies_.SetL(L_0th);

        // initialize the true p and L (i.e. including the perturbations)
        VectorArray<N> p = p_of_v();
        QuaternionArray<N> L = L_of_w();
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

    // get all current bodies
    const BodyArray<N>& GetCurrentBodies() const {
        return bodies_;
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

    heap_array<double, N> lambda0_;

    VectorArray<N> x_cached_Phi_;
    VectorArray<N> orient_cached_I_;


    VectorArray<N> p_of_v()
    {
        // get the momenta from the velocities
        Compute_Phi_grad_Phi();
        Compute_Theta_grad_Theta();

        VectorArray<N> p;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& M = body.GetM();
            const Vector& v = body.Getv();

            p.row(i) = M*v
                        + PHYS_inv_c2 *(M/2 *v.squaredNorm()*v - 3*M *Phi_[i]*v - 4*M *Theta_[i]);
        }

        return p;
    }


    VectorArray<N> v_of_p()
    {
        // get the velocities from the momenta
        Compute_Phi_grad_Phi();
        Compute_Theta_grad_Theta();

        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& M = body.GetM();
            const Vector& p = body.Getp();

            v.row(i) = p/M
                        + PHYS_inv_c2 *(-1.0/(2*M*M*M) *p.squaredNorm()*p + 3.0/M *Phi_[i]*p + 4.0 *Theta_[i]);
        }

        return v;
    }


    QuaternionArray<N> L_of_w()
    {
        // get the angular momenta from the angular velocities
        QuaternionArray<N> L;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& Ixy = body.GetIxy();
            const double& Iz = body.GetIz();
            const Quaternion& w = body.Getw();
            const Quaternion& q = body.Getq();

            double Re = ((q.conjugate() *w).mult_k()).scalar();

            L.row(i) = 4*Ixy*w - 4*(Iz - Ixy) *Re *q.mult_k();
        }

        return L;
    }


    QuaternionArray<N> w_of_L()
    {
        // get the angular velocities from the angular momenta
        QuaternionArray<N> w;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& Ixy = body.GetIxy();
            const double& Iz = body.GetIz();
            const Quaternion& L = body.GetL();
            const Quaternion& q = body.Getq();

            double Re = ((q.conjugate() *L).mult_k()).scalar();

            w.row(i) = 1.0/(4*Ixy)*L + 1.0/(4*Ixy*Ixy)*(Iz - Ixy) *Re *q.mult_k();
        }

        return w;
    }


    VectorArray<N> dH_kin_dp()
    {
        // derivative of the kinetic Hamiltonian wrt. momentum
        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const Vector& p = body.Getp();
            const double& M = body.GetM();

            v.row(i) = p/M
                        + PHYS_inv_c2 *(-1.0/(2*M*M*M) *p.squaredNorm()*p);
        }

        return v;
    }


    VectorArray<N> dH_pot_dx()
    {
        // derivative of the potential Hamiltonian wrt. position
        Compute_Phi_grad_Phi();

        VectorArray<N> NegF;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& M = body.GetM();

            NegF.row(i) = M*grad_Phi_[i]
                            + PHYS_inv_c2 *(M *Phi_[i]*grad_Phi_[i]);
        }

        for (unsigned i = 0; i < N; i++) {
            const Body& body_x = bodies_[i];
            const Vector& x = body_x.Getx();
            const double& Mx = body_x.GetM();
            const double& ax = body_x.Geta();
            const double& J2x = body_x.GetJ2();
            const Vector epx = body_x.GetPoleAxis();

            for (unsigned j = 0; j < i; j++) {
                const Body& body_y = bodies_[j];
                const Vector& y = body_y.Getx();
                const double& My = body_y.GetM();
                const double& ay = body_y.Geta();
                const double& J2y = body_y.GetJ2();
                const Vector epy = body_y.GetPoleAxis();

                Vector r_vec = x - y;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;
                double r2_inv = 1/(r*r);
                double epxer = epx.dot(e_r);
                double epyer = epy.dot(e_r);

                Vector negf = PHYS_G*Mx*My*r2_inv *(
                    -3*J2x*(ax*ax*r2_inv)*((2.5*epxer*epxer - 0.5)*e_r - epxer*epx)
                    -3*J2y*(ay*ay*r2_inv)*((2.5*epyer*epyer - 0.5)*e_r - epyer*epy));

                NegF.row(i) += negf;
                NegF.row(j) -= negf;
            }
        }

        return NegF;
    }


    QuaternionArray<N> dH_kin_dL()
    {
        // derivative of the kinetic Hamiltonian wrt. angular momentum
        QuaternionArray<N> w;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const Quaternion& L = body.GetL();
            const double& Ixy = body.GetIxy();

            w.row(i) = L /(4*Ixy);
        }

        return w;
    }


    QuaternionArray<N> dH_pot_dq()
    {
        // derivative of the potential Hamiltonian wrt. orientation
        QuaternionArray<N> NegTorque;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const Quaternion& q = body.Getq();

            NegTorque.row(i) = 2*lambda0_[i] *q;
        }

        for (unsigned i = 0; i < N; i++) {
            const Body& body_x = bodies_[i];
            const Vector& x = body_x.Getx();
            const double& Mx = body_x.GetM();
            const double& ax = body_x.Geta();
            const double& J2x = body_x.GetJ2();
            const Vector epx = body_x.GetPoleAxis();
            const Quaternion& qx = body_x.Getq();

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}

                const Body& body_y = bodies_[j];
                const Vector& y = body_y.Getx();
                const double& My = body_y.GetM();

                Vector r_vec = x - y;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;
                double r3_inv = 1/(r*r*r);
                double epxer = epx.dot(e_r);

                Quaternion negtorque = -6*PHYS_G*Mx*My*r3_inv *J2x*(ax*ax) *epxer *(e_r*qx).mult_k();

                NegTorque.row(i) += negtorque;
            }
        }

        return NegTorque;
    }


    VectorArray<N> dH_pert_dx()
    {
        // derivative of the perturbation Hamiltonian wrt. position
        Compute_Phi_grad_Phi();
        Compute_Theta_grad_Theta();
        Compute_grad_Lambda();

        VectorArray<N> NegF;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& M = body.GetM();
            const Vector& p = body.Getp();

            NegF.row(i) = PHYS_inv_c2 *(4 *grad_Theta_[i]*p + 3.0/(2*M) *grad_Phi_[i]*p.squaredNorm() + M/2 *grad_Lambda_[i]);
        }

        return NegF;
    }


    VectorArray<N> dH_pert_dp()
    {
        // derivative of the perturbation Hamiltonian wrt. momentum
        Compute_Phi_grad_Phi();
        Compute_Theta_grad_Theta();

        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& M = body.GetM();
            const Vector& p = body.Getp();

            v.row(i) = PHYS_inv_c2 *(4.0 *Theta_[i] + 3.0/M *Phi_[i]*p);
        }

        return v;
    }


    QuaternionArray<N> dH_pert_dq()
    {
        // derivative of the perturbation Hamiltonian wrt. orientation

        QuaternionArray<N> NegTorque;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& Ixy = body.GetIxy();
            const double& Iz = body.GetIz();
            const Quaternion& L = body.GetL();
            const Quaternion& q = body.Getq();

            double lambda = 1.0/(8*Ixy)* L.squaredNorm();
            double Re = ((q.conjugate() *L).mult_k()).scalar();

            NegTorque.row(i) = -1.0/(4*Ixy*Ixy)*(Iz - Ixy) *Re *L.mult_k()
                            + 2*(lambda - lambda0_[i])*q;
        }

        return NegTorque;
    }


    QuaternionArray<N> dH_pert_dL()
    {
        // derivative of the perturbation Hamiltonian wrt. angular momentum
        QuaternionArray<N> w;
        for (unsigned i = 0; i < N; i++) {
            const Body& body = bodies_[i];
            const double& Ixy = body.GetIxy();
            const double& Iz = body.GetIz();
            const Quaternion& L = body.GetL();
            const Quaternion& q = body.Getq();

            double Re = ((q.conjugate() *L).mult_k()).scalar();

            w.row(i) = 1.0/(4*Ixy*Ixy)*(Iz - Ixy) *Re *q.mult_k();
        }

        return w;
    }


    void Compute_Phi_grad_Phi() {
        // compute the leading order metric perturbation Phi (aka. the Newtonian gravitational potential)
        // as well as its gradient
        if (bodies_.Getx() == x_cached_Phi_) {return;}

        for (unsigned i = 0; i < N; i++) {
            Phi_[i] = 0.0;
            grad_Phi_[i] = Vector::Zero();

            const Vector& x0 = bodies_[i].Getx();

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}

                const Body& body = bodies_[j];
                const double& M = body.GetM();
                const Vector& x = body.Getx();

                Vector r_vec = x0 - x;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;

                Phi_[i] += -PHYS_G*M /r;
                grad_Phi_[i] += PHYS_G*M /(r*r) *e_r;
            }
        }

        x_cached_Phi_ = bodies_.Getx();
    }


    void Compute_Theta_grad_Theta() {
        // compute the metric perturbation Theta
        // as well as its gradient

        for (unsigned i = 0; i < N; i++) {
            Theta_[i] = Vector::Zero();
            grad_Theta_[i] = Matrix::Zero();

            const Vector& x0 = bodies_[i].Getx();

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}

                const Body& body = bodies_[j];
                const double& M = body.GetM();
                const Vector& x = body.Getx();
                const Vector& p = body.Getp();

                Vector v = p/M;

                Vector r_vec = x0 - x;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;

                Theta_[i] += PHYS_G*M /r *v;
                grad_Theta_[i] += -PHYS_G*M /(r*r) *e_r *v.transpose();
            }
        }
    }


    void Compute_grad_Lambda() {
        // compute the gradient of metric perturbation Lambda
        Compute_Phi_grad_Phi();

        for (unsigned i = 0; i < N; i++) {
            grad_Lambda_[i] = Vector::Zero();

            const Vector& x0 = bodies_[i].Getx();

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}

                const Body& body = bodies_[j];
                const double& M = body.GetM();
                const Vector& x = body.Getx();
                const Vector& p = body.Getp();

                Vector r_vec = x0 - x;
                double r = r_vec.norm();
                Vector e_r = r_vec /r;

                Vector v = p/M;
                Vector a = -grad_Phi_[j];

                Vector incr_1 = 4*PHYS_G*M *v.squaredNorm() /(r*r) * e_r;
                Vector incr_2 = PHYS_G*M* e_r.dot(v) /(r*r) *(2*v - 3*e_r.dot(v) *e_r);
                Vector incr_3 = -4*PHYS_G*M *Phi_[j] /(r*r) * e_r;
                Vector incr_4 = PHYS_G*M /r *(a - (e_r.dot(a)) *e_r);

                grad_Lambda_[i] += incr_1 + incr_2 + incr_3 + incr_4;
            }
        }
    }

};









#endif // PHYSICS_H_INCLUDED

