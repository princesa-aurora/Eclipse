#include <iostream>
#include <functional>
#include <vector>
#include <chrono>
#include <cstdio>

#include <include/utils.h>


template<unsigned N>
VectorArray<N> dH_kin_dp(const BodyArray<N> &bodies)
{
    // derivative of the kinetic Hamiltonian wrt. momentum
    VectorArray<N> v;
    for (unsigned i = 0; i < N; i++) {
        const Body& body = bodies[i];
        const Vector& x = body.Getx();
        const Vector& p = body.Getp();
        const double& M = body.GetM();
        v.row(i) = body.Getv();

        v.row(i) *= (1.0 - p.squaredNorm()/(2*M*M*PHYS_c*PHYS_c));
    }

    return v;
};


template<unsigned N>
VectorArray<N> dH_pot_dx(const BodyArray<N> &bodies)
{
    // derivative of the potetntial Hamiltonian wrt. position
    const Vector& x_sun = bodies[0].Getx();
    const double& M_sun = bodies[0].GetM();

    VectorArray<N> NegF = VectorArray<N>::Zero();
    Vector negf;
    Vector r_vec;
    double r;
    Vector e_r;
    double r2_inv;
    double ewxer;
    double ewyer;

    for (unsigned i = 0; i < N; i++) {
        const Vector& x = bodies[i].Getx();
        const double& Mx = bodies[i].GetM();
        const double& ax = bodies[i].Geta();
        const double& J2x = bodies[i].GetJ2();
        const Vector ewx = bodies[i].GetAxis();

        for (unsigned j = 0; j < i; j++) {
            const Vector& y = bodies[j].Getx();
            const double& My = bodies[j].GetM();
            const double& ay = bodies[j].Geta();
            const double& J2y = bodies[j].GetJ2();
            const Vector ewy = bodies[j].GetAxis();

            r_vec = y - x;
            r = r_vec.norm();
            e_r = r_vec /r;
            r2_inv = 1/(r*r);
            ewxer = ewx.dot(e_r);
            ewyer = ewy.dot(e_r);

            negf = -PHYS_G*Mx*My*r2_inv *(e_r
                +3*J2x*(ax*ax*r2_inv)*((2.5*ewxer*ewxer - 0.5)*e_r - ewxer*ewx)
                +3*J2y*(ay*ay*r2_inv)*((2.5*ewyer*ewyer - 0.5)*e_r - ewyer*ewy));

            NegF.row(i) += negf;
            NegF.row(j) -= negf;
        }
    }

    for (unsigned i = 1; i < N; i++) {
        const Body& body = bodies[i];
        const Vector& x = body.Getx();
        const Vector& p = body.Getp();
        const double& M = body.GetM();
        //NegF.row(i) *= (1.0 + PHYS_G*M_sun/((x-x_sun).norm()*PHYS_c*PHYS_c));
    }
    return NegF;
};


template<unsigned N>
VectorArray<N> dH_kin_dL(const BodyArray<N> &bodies)
{
    // derivative of the kinetic Hamiltonian wrt. angular momentum
    VectorArray<N> w;
    for (unsigned i = 0; i < N; i++) {
        w.row(i) = bodies[i].Getw();
    }
    return w;
};


template<unsigned N>
VectorArray<N> dH_pot_dorient(const BodyArray<N> &bodies)
{
    // derivative of the potential Hamiltonian wrt. orientation
    VectorArray<N> NegTorque = VectorArray<N>::Zero();
    Vector negtorque;
    Vector r_vec;
    double r;
    Vector e_r;
    double r3_inv;
    double ewer;

    for (unsigned i = 0; i < N; i++) {
        const Vector& x = bodies[i].Getx();
        const double& alpha = bodies[i].Getorient()(1);
        const double& delta = bodies[i].Getorient()(2);
        const Vector& L = bodies[i].GetL();
        const double& Mx = bodies[i].GetM();
        const double& Ixy = bodies[i].GetIxy();
        const double& Iz = bodies[i].GetIz();
        const double& a = bodies[i].Geta();
        const double& J2 = bodies[i].GetJ2();
        const Vector ew = bodies[i].GetAxis();

        const Vector ealpha(-sin(alpha), cos(alpha), 0.0);
        const Vector edelta(-sin(delta)*cos(alpha), -sin(delta)*sin(alpha), cos(delta));

        Matrix I_inv_prime;
        I_inv_prime << -2.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 0.0,
                        -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 2.0/Ixy*sin(delta)/cos(delta)/cos(delta)/cos(delta), 0.0,
                        0.0, 0.0, 0.0;

        for (unsigned j = 0; j < N; j++) {
            if (j == i) {continue;}
            const Vector& y = bodies[j].Getx();
            const double& My = bodies[j].GetM();

            r_vec = y - x;
            r = r_vec.norm();
            e_r = r_vec /r;
            r3_inv = 1/(r*r*r);
            ewer = ew.dot(e_r);

            negtorque(0) = 0.0;
            negtorque(1) = -3*PHYS_G*Mx*My*r3_inv*a*a*J2*ewer*(ealpha.dot(e_r))*cos(delta);
            negtorque(2) = 1.0/2*(I_inv_prime *L).dot(L)
                        - 3*PHYS_G*Mx*My*r3_inv*a*a*J2*ewer*(edelta.dot(e_r));

            NegTorque.row(i) += negtorque;
        }
    }
    return NegTorque;
};


template<unsigned N>
VectorArray<N> dH_pert_dx(const BodyArray<N> &bodies)
{
    // derivative of the perturbation Hamiltonian wrt. position
    return VectorArray<N>::Zero();
};

template<unsigned N>
VectorArray<N> dH_pert_dp(const BodyArray<N> &bodies)
{
    // derivative of the perturbation Hamiltonian wrt. momentum
    return VectorArray<N>::Zero();
};

