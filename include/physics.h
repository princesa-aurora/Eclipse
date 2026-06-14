#include <iostream>
#include <functional>
#include <vector>
#include <chrono>
#include <cstdio>

#include <include/utils.h>


template<unsigned N>
VectorArray<N> p_of_v(const BodyArray<N> &bodies)
{
    // get the momenta from the velocities
    VectorArray<N> p;
    const double& M_sun = bodies[0].GetM();
    const Vector x_sun = bodies[0].Getx();
    for (unsigned i = 0; i < N; i++) {
        const Vector& v = bodies[i].Getv();
        const Vector& x = bodies[i].Getx();
        const double& M = bodies[i].GetM();
        p.row(i) = M*v * (1 + v.squaredNorm()/(2*PHYS_c*PHYS_c));
        if (i != 0) {
            p.row(i) *= (1 + 3*PHYS_G*M_sun/(PHYS_c*PHYS_c*(x-x_sun).norm()));
        }
    }

    return p;
}


template<unsigned N>
VectorArray<N> v_of_p(const BodyArray<N> &bodies)
{
    // get the velocities from the momenta
    VectorArray<N> v;
    const double& M_sun = bodies[0].GetM();
    const Vector x_sun = bodies[0].Getx();
    for (unsigned i = 0; i < N; i++) {
        const Vector& p = bodies[i].Getp();
        const Vector& x = bodies[i].Getx();
        const double& M = bodies[i].GetM();
        v.row(i) = p/M * (1 - p.squaredNorm()/(2*PHYS_c*PHYS_c*M*M));
        if (i != 0) {
            v.row(i) *= (1 - 3*PHYS_G*M_sun/(PHYS_c*PHYS_c*(x-x_sun).norm()));
        }
    }

    return v;
}


template<unsigned N>
VectorArray<N> L_of_w(const BodyArray<N> &bodies)
{
    // get the angular momenta from the angular velocities
    VectorArray<N> L;
    for (unsigned i = 0; i < N; i++) {
        const Vector& w = bodies[i].Getw();
        const double& delta = bodies[i].Getorient()(2);
        const double& Ixy = bodies[i].GetIxy();
        const double& Iz = bodies[i].GetIz();
        Matrix I;
        I << Iz, -Iz*(1.0-sin(delta)), 0.0,
             -Iz*(1.0-sin(delta)), Ixy*cos(delta)*cos(delta) + Iz*(1.0-sin(delta))*(1.0-sin(delta)), 0.0,
             0.0, 0.0,Ixy;
        L.row(i) = I *w;
    }

    return L;
}


template<unsigned N>
VectorArray<N> w_of_L(const BodyArray<N> &bodies)
{
    // get the angular velocities from the angular momenta
    VectorArray<N> w;
    for (unsigned i = 0; i < N; i++) {
        const Vector& L = bodies[i].GetL();
        const double& delta = bodies[i].Getorient()(2);
        const double& Ixy = bodies[i].GetIxy();
        const double& Iz = bodies[i].GetIz();
        Matrix I_inv;
        I_inv << 1.0/Iz + 1.0/Ixy*cos(delta)*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 1.0/Ixy/(1.0+sin(delta)), 0.0,
                 1.0/Ixy/(1.0+sin(delta)), 1.0/Ixy/cos(delta)/cos(delta), 0.0,
                 0.0, 0.0, 1.0/Ixy;
        w.row(i) = I_inv *L;
    }

    return w;
}


template<unsigned N>
VectorArray<N> dH_kin_dp(const BodyArray<N> &bodies)
{
    // derivative of the kinetic Hamiltonian wrt. momentum
    VectorArray<N> v;
    for (unsigned i = 0; i < N; i++) {
        const Vector& p = bodies[i].Getp();
        const double& M = bodies[i].GetM();
        v.row(i) = p/M * (1.0 - p.squaredNorm()/(2*M*M*PHYS_c*PHYS_c));
    }

    return v;
};


template<unsigned N>
VectorArray<N> dH_pot_dx(const BodyArray<N> &bodies)
{
    // derivative of the potetntial Hamiltonian wrt. position
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

            r_vec = x - y;
            r = r_vec.norm();
            e_r = r_vec /r;
            r2_inv = 1/(r*r);
            ewxer = ewx.dot(e_r);
            ewyer = ewy.dot(e_r);

            negf = PHYS_G*Mx*My*r2_inv *(e_r
                -3*J2x*(ax*ax*r2_inv)*((2.5*ewxer*ewxer - 0.5)*e_r - ewxer*ewx)
                -3*J2y*(ay*ay*r2_inv)*((2.5*ewyer*ewyer - 0.5)*e_r - ewyer*ewy));

            NegF.row(i) += negf;
            NegF.row(j) -= negf;

            if (j == 0) {
                NegF.row(i) += -PHYS_G*PHYS_G*My*My*Mx/(PHYS_c*PHYS_c*r*r*r) *e_r;
            }
        }
    }

    return NegF;
};


template<unsigned N>
VectorArray<N> dH_kin_dL(const BodyArray<N> &bodies)
{
    // derivative of the kinetic Hamiltonian wrt. angular momentum
    return w_of_L(bodies);
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
    const Vector& x_sun = bodies[0].Getx();
    const double& M_sun = bodies[0].GetM();

    VectorArray<N> NegF;
    Vector r_vec;
    double r;
    Vector e_r;

    NegF.row(0) = Vector::Zero();

    for (unsigned i = 1; i < N; i++) {
        const Vector& x = bodies[i].Getx();
        const Vector& p = bodies[i].Getp();
        const double& Mx = bodies[i].GetM();

        r_vec = x - x_sun;
        r = r_vec.norm();
        e_r = r_vec /r;

        NegF.row(i) = 3*PHYS_G*M_sun/(2*Mx*PHYS_c*PHYS_c) *p.squaredNorm() /(r*r) *e_r;
    }

    return NegF;
};

template<unsigned N>
VectorArray<N> dH_pert_dp(const BodyArray<N> &bodies)
{
    // derivative of the perturbation Hamiltonian wrt. momentum
    const Vector& x_sun = bodies[0].Getx();
    const double& M_sun = bodies[0].GetM();

    VectorArray<N> v;
    Vector r_vec;
    double r;
    Vector e_r;

    v.row(0) = Vector::Zero();

    for (unsigned i = 1; i < N; i++) {
        const Vector& x = bodies[i].Getx();
        const Vector& p = bodies[i].Getp();
        const double& Mx = bodies[i].GetM();

        r_vec = x - x_sun;
        r = r_vec.norm();
        e_r = r_vec /r;

        v.row(i) = -3*PHYS_G*M_sun/(Mx*PHYS_c*PHYS_c) *p /r;
    }

    return v;
};

