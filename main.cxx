#include <iostream>
#include <functional>
#include <vector>

#include <include/utils.h>
#include <include/constants.h>
#include <include/initial_conditions.h>
#include <include/solver.h>

//#define WITHOUT_NUMPY
//#include <include/matplotlibcpp.h>
//namespace plt = matplotlibcpp;

const int N = 4; // number of bodies in the system

using VectorArray = Array<Vector, N>;
using ScalarArray = Array<double, N>;

double t0 = t0_mjdtdb;
VectorArray x0 = std::array{x0_sun, x0_earth, x0_moon, x0_jupiter};
VectorArray p0 = std::array{p0_sun, p0_earth, p0_moon, p0_jupiter};
ScalarArray M = std::array{M_sun, M_earth, M_moon, M_jupiter};
ScalarArray GM = std::array{GM_sun, GM_earth, GM_moon, GM_jupiter};

double dt = 60.0;
double T = 200.0 * 365.25 * 24 * 60 * 60;
int K = T/dt;



int main() {

    std::function<VectorArray(const VectorArray&, const VectorArray&)> x_dot = [&](const VectorArray &x, const VectorArray &p)
    {
        return p/M;
    };

    std::function<VectorArray(const VectorArray&, const VectorArray&)> p_dot = [&](const VectorArray &x, const VectorArray &p)
    {
        VectorArray F;
        Vector dist;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (j == i) {continue;}

                dist = x.at(j) - x.at(i);
                F[i] = F.at(i) + GM.at(j)*M.at(i)*dist/pow(dist*dist, 3.0/2);
            }
        }

        return F;
    };

    Forest_Ruth solver(x_dot, p_dot, x0, p0, t0, false);

    double t = 0.0;
    VectorArray pos;
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


    std::cout << "Starting integration with dt = " << dt << " seconds, for a total time of " << T << " seconds (" << K << " steps)." << "\n" << std::endl;

    for (int i = 0; i < K; i++) {
        solver.MakeStep(dt);

        t = solver.GetCurrentTime();
        pos = solver.GetCurrentPosition();

        lunar_eclipse = eclipsed(pos.at(2), R_moon, pos.at(1), R_earth, pos.at(0), R_sun);
        solar_eclipse = eclipsed(pos.at(0), R_sun, pos.at(2), R_moon, pos.at(1), R_earth);

        if (lunar_eclipse && !lunar_eclipse_active) {
            std::cout << "detected lunar eclipse start at time " << mjdtdb_to_utc_date(t/86400.0) << " seconds." << std::endl;
        }
        if (!lunar_eclipse && lunar_eclipse_active) {
            std::cout << "detected lunar eclipse maximum at time " << mjdtdb_to_utc_date(lunar_eclipse_max_t/86400.0) << " with occultation " << int(lunar_eclipse_max_ratio*1000)/10.0 << "%." << std::endl;
            std::cout << "detected lunar eclipse end at time " << mjdtdb_to_utc_date(t/86400.0) << " seconds." << "\n" << std::endl;
            lunar_eclipse_max_ratio = 0.0;
        }
        lunar_eclipse_active = lunar_eclipse;

        if (lunar_eclipse_active) {
            lunar_eclipse_ratio = eclipse_ratio(pos.at(2), R_moon, pos.at(1), R_earth, pos.at(0), R_sun);

            if (lunar_eclipse_ratio > lunar_eclipse_max_ratio) {
                lunar_eclipse_max_ratio = lunar_eclipse_ratio;
                lunar_eclipse_max_t = t;
            }
        }

        if (solar_eclipse && !solar_eclipse_active) {
            std::cout << "detected solar eclipse start at time " << mjdtdb_to_utc_date(t/86400.0) << " seconds." << std::endl;
        }
        if (!solar_eclipse && solar_eclipse_active) {
            std::cout << "detected solar eclipse maximum at time " << mjdtdb_to_utc_date(solar_eclipse_max_t/86400.0) << " with occultation " << int(solar_eclipse_max_ratio*1000)/10.0 << "%." << std::endl;
            std::cout << "detected solar eclipse end at time " << mjdtdb_to_utc_date(t/86400.0) << " seconds." << "\n" << std::endl;
            solar_eclipse_max_ratio = 0.0;
        }
        solar_eclipse_active = solar_eclipse;

        if (solar_eclipse_active) {
            solar_eclipse_ratio = eclipse_ratio(pos.at(0), R_sun, pos.at(2), R_moon, pos.at(1), R_earth);

            if (solar_eclipse_ratio > solar_eclipse_max_ratio) {
                solar_eclipse_max_ratio = solar_eclipse_ratio;
                solar_eclipse_max_t = t;
            }
        }
    }

    std::cout << "Integration complete." << std::endl;
}



