#include <iostream>
#include <functional>
#include <vector>

#include <include/utils.h>
#include <include/initial_conditions.h>
#include <include/solver.h>

//#define WITHOUT_NUMPY
//#include <include/matplotlibcpp.h>
//namespace plt = matplotlibcpp;

const size_t N = 10; // number of bodies in the system

double t0 = 0.0; // seconds since J2000.0 epoch
BodyArray<N> initial_bodies(std::array{Sun, Earth, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune});

double dt = 60.0;
double T = 50.0 * 365.25 * 24 * 60 * 60;
size_t K = T/dt;



int main() {

    std::function<VectorArray<N>(const BodyArray<N>&)> x_dot = [&](const BodyArray<N> &bodies)
    {
        // position time derivative
        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            v[i] = bodies.Getv(i);
        }
        return v;
    };

    std::function<VectorArray<N>(const BodyArray<N>&)> p_dot = [&](const BodyArray<N> &bodies)
    {
        // momentum time derivative
        VectorArray<N> F;
        Vector f;
        Vector r_vec;
        double r;
        Vector e_r;
        Matrix proj;

        for (unsigned i = 0; i < N; i++) {
            const Vector& x = bodies.Getx(i);
            const double& Mx = bodies.GetM(i);
            const double& Rx = bodies.GetR(i);
            const double& J2x = bodies.GetJ2(i);
            const Vector ewx = bodies.GetAxis(i);

            for (unsigned j = 0; j < i; j++) {
                const Vector& y = bodies.Getx(j);
                const double& My = bodies.GetM(j);
                const double& Ry = bodies.GetR(j);
                const double& J2y = bodies.GetJ2(j);
                const Vector ewy = bodies.GetAxis(j);

                r_vec = y - x;
                r = r_vec.Abs();
                e_r = r_vec /r;
                proj = Identity - TensorProduct(e_r, e_r);

                f = PHYS_G*Mx*My/r/r *(e_r
                    +3*J2x*(Rx*Rx/r/r)*(P2(ewx*e_r)*e_r - (ewx*e_r)*proj*ewx)
                    +3*J2y*(Ry*Ry/r/r)*(P2(ewy*e_r)*e_r - (ewy*e_r)*proj*ewy));

                F[i] += f;
                F[j] -= f;
            }
        }

        return F;
    };


    std::function<VectorArray<N>(const BodyArray<N>&)> orient_dot = [&](const BodyArray<N> &bodies)
    {
        // rotation angle time derivative
        VectorArray<N> w;
        for (unsigned i = 0; i < N; i++) {
            w[i] = bodies.Getw(i);
        }
        return w;
    };

    std::function<VectorArray<N>(const BodyArray<N>&)> L_dot = [&](const BodyArray<N> &bodies)
    {
        // angular momentum time derivative
        VectorArray<N> Torque;
        Vector torque;
        Vector r_vec;
        double r;
        Vector e_r;

        for (unsigned i = 0; i < N; i++) {
            const Vector& x = bodies.Getx(i);
            const double& alpha = bodies.Getorient(i)[1];
            const double& delta = bodies.Getorient(i)[2];
            const Vector& L = bodies.GetL(i);
            const double& Mx = bodies.GetM(i);
            const double& Ixy = bodies.GetIxy(i);
            const double& Iz = bodies.GetIz(i);
            const double& R = bodies.GetR(i);
            const double& J2 = bodies.GetJ2(i);
            const Vector ew = bodies.GetAxis(i);

            const Vector ealpha{-sin(alpha), cos(alpha), 0.0};
            const Vector edelta{-sin(delta)*cos(alpha), -sin(delta)*sin(alpha), cos(delta)};

            const Matrix I_inv_prime{{-2.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 0.0},
                                     {-1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 2.0/Ixy*sin(delta)/cos(delta)/cos(delta)/cos(delta), 0.0},
                                     {0.0, 0.0, 0.0}};

            for (unsigned j = 0; j < i; j++) {
                const Vector& y = bodies.Getx(j);
                const double& My = bodies.GetM(j);

                r_vec = y - x;
                r = r_vec.Abs();
                e_r = r_vec /r;

                torque[0] = 0.0;
                torque[1] = 3*PHYS_G*Mx*My/r/r/r*R*R*J2*(ew*e_r)*(ealpha*e_r)*cos(delta);
                torque[2] = -1.0/2*L* (I_inv_prime *L)
                            + 3*PHYS_G*Mx*My/r/r/r*R*R*J2*(ew*e_r)*(edelta*e_r);

                Torque[i] += torque;
                Torque[j] -= torque;
            }
        }

        return Torque;
    };



    Forest_Ruth<N> solver(x_dot, p_dot, orient_dot, L_dot, initial_bodies, t0);

    double t = 0.0;
    Vector x_sun;
    Vector x_earth;
    Vector x_moon;
    bool lunar_eclipse = false;
    bool solar_eclipse = false;
    bool lunar_eclipse_active = false;
    bool solar_eclipse_active = false;
    double lunar_eclipse_ratio = 0.0;
    double solar_eclipse_ratio = 0.0;
    double lunar_eclipse_max_ratio = 0.0;
    double solar_eclipse_max_ratio = 0.0;
    double lunar_eclipse_max_t = 0.0;
    double solar_eclipse_max_t = 0.0;

    double RA_Erlangen = 11.0*M_PI/180;
    double Dec_Erlangen = 49.6*M_PI/180;
    double zenith;
    bool day = true;
    bool day_active = true;


    std::cout << "Starting integration with dt = " << dt << " seconds for a total time of " << T << " seconds (" << T/60/60/24/365.25 << " years) (" << K << " steps)." << "\n" << std::endl;

    for (size_t i = 0; i < K; i++) {
        solver.MakeStep(dt);

        t = solver.GetCurrentTime();
        x_sun = solver.GetCurrentBody(0).Getx();
        x_earth = solver.GetCurrentBody(1).Getx();
        x_moon = solver.GetCurrentBody(2).Getx();

        lunar_eclipse = eclipsed(x_moon, R_moon, x_earth, R_earth, x_sun, R_sun);
        solar_eclipse = eclipsed(x_sun, R_sun, x_moon, R_moon, x_earth, R_earth);

        if (lunar_eclipse && !lunar_eclipse_active) {
            std::cout << "detected lunar eclipse start at time " << j2000_to_utc_date(t) << "."<< std::endl;
        }
        if (!lunar_eclipse && lunar_eclipse_active) {
            std::cout << "detected lunar eclipse maximum at time " << j2000_to_utc_date(lunar_eclipse_max_t) << " with occultation " << int(lunar_eclipse_max_ratio*1000)/10.0 << "%." << std::endl;
            std::cout << "detected lunar eclipse end at time " << j2000_to_utc_date(t) << "." << "\n" << std::endl;
            lunar_eclipse_max_ratio = 0.0;
        }
        lunar_eclipse_active = lunar_eclipse;

        if (lunar_eclipse_active) {
            lunar_eclipse_ratio = eclipse_ratio(x_moon, R_moon, x_earth, R_earth, x_sun, R_sun);

            if (lunar_eclipse_ratio > lunar_eclipse_max_ratio) {
                lunar_eclipse_max_ratio = lunar_eclipse_ratio;
                lunar_eclipse_max_t = t;
            }
        }

        if (solar_eclipse && !solar_eclipse_active) {
            std::cout << "detected solar eclipse start at time " << j2000_to_utc_date(t) << "." << std::endl;
        }
        if (!solar_eclipse && solar_eclipse_active) {
            std::cout << "detected solar eclipse maximum at time " << j2000_to_utc_date(solar_eclipse_max_t) << " with occultation " << int(solar_eclipse_max_ratio*1000)/10.0 << "%." << std::endl;
            std::cout << "detected solar eclipse end at time " << j2000_to_utc_date(t) << "." << "\n" << std::endl;
            solar_eclipse_max_ratio = 0.0;
        }
        solar_eclipse_active = solar_eclipse;

        if (solar_eclipse_active) {
            solar_eclipse_ratio = eclipse_ratio(x_sun, R_sun, x_moon, R_moon, x_earth, R_earth);

            if (solar_eclipse_ratio > solar_eclipse_max_ratio) {
                solar_eclipse_max_ratio = solar_eclipse_ratio;
                solar_eclipse_max_t = t;
            }
        }

        zenith = sun_zenith(RA_Erlangen, Dec_Erlangen, solver.GetCurrentBody(1), solver.GetCurrentBody(0));
        day = (zenith > 0.0);

        if (day && !day_active) {
            std::cout << "Erlangen, Germany: sunrise detected at time " << j2000_to_utc_date(t) << ".\n" << std::endl;
        }
        if (!day && day_active) {
            std::cout << "Erlangen, Germany: sunset detected at time " << j2000_to_utc_date(t) << ".\n" << std::endl;
        }
        day_active = day;

    }

    std::cout << "Integration complete." << std::endl;
}



