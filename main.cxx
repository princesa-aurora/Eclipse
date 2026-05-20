#include <iostream>
#include <functional>
#include <vector>
#include <chrono>
#include <filesystem>
#include <format>

#include <include/utils.h>
#include <include/initial_conditions.h>
#include <include/solver.h>

namespace fs = std::filesystem;


const size_t N = 10; // number of bodies in the system

double t0 = 0.0; // seconds since J2000.0 epoch
BodyArray<N> initial_bodies(std::array{Sun, Earth, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune});

double dt = 60.0;
double T = 50.0 * 365.25 * 24 * 60 * 60;
size_t K = T/dt;



int main() {
    // name and create output folder using a timestamp
    std::string timestamp = std::format("{:%d.%m.%Y_%H:%M:%S}",
                                        std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::system_clock::now()));
    fs::path folder = "/home/aurora/eclipse_data" / fs::path(timestamp);
    fs::create_directory(folder);


    std::function<VectorArray<N>(const BodyArray<N>&)> x_dot = [&](const BodyArray<N> &bodies)
    {
        // position time derivative
        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            v.row(i) = bodies.Getv(i);
        }
        return v;
    };

    std::function<VectorArray<N>(const BodyArray<N>&)> p_dot = [&](const BodyArray<N> &bodies)
    {
        // momentum time derivative

        //std::chrono::_V2::system_clock::time_point start = std::chrono::high_resolution_clock::now();

        VectorArray<N> F = VectorArray<N>::Zero();
        Vector f;
        Vector r_vec;
        double r;
        Vector e_r;
        double r2_inv;
        double ewxer;
        double ewyer;

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
                r = r_vec.norm();
                e_r = r_vec /r;
                r2_inv = 1/(r*r);
                ewxer = ewx.dot(e_r);
                ewyer = ewy.dot(e_r);

                f = PHYS_G*Mx*My*r2_inv *(e_r
                    +3*J2x*(Rx*Rx*r2_inv)*((2.5*ewxer*ewxer - 0.5)*e_r - ewxer*ewx)
                    +3*J2y*(Ry*Ry*r2_inv)*((2.5*ewyer*ewyer - 0.5)*e_r - ewyer*ewy));

                F.row(i) += f;
                F.row(j) -= f;
            }
        }

        //std::chrono::_V2::system_clock::time_point stop = std::chrono::high_resolution_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;

        return F;
    };


    std::function<VectorArray<N>(const BodyArray<N>&)> orient_dot = [&](const BodyArray<N> &bodies)
    {
        // rotation angle time derivative
        VectorArray<N> w;
        for (unsigned i = 0; i < N; i++) {
            w.row(i) = bodies.Getw(i);
        }
        return w;
    };

    std::function<VectorArray<N>(const BodyArray<N>&)> L_dot = [&](const BodyArray<N> &bodies)
    {
        // angular momentum time derivative

        //std::chrono::_V2::system_clock::time_point start = std::chrono::high_resolution_clock::now();

        VectorArray<N> Torque = VectorArray<N>::Zero();
        Vector torque;
        Vector r_vec;
        double r;
        Vector e_r;
        double r3_inv;
        double ewer;

        for (unsigned i = 0; i < N; i++) {
            const Vector& x = bodies.Getx(i);
            const double& alpha = bodies.Getorient(i)(1);
            const double& delta = bodies.Getorient(i)(2);
            const Vector& L = bodies.GetL(i);
            const double& Mx = bodies.GetM(i);
            const double& Ixy = bodies.GetIxy(i);
            const double& Iz = bodies.GetIz(i);
            const double& R = bodies.GetR(i);
            const double& J2 = bodies.GetJ2(i);
            const Vector ew = bodies.GetAxis(i);

            const Vector ealpha(-sin(alpha), cos(alpha), 0.0);
            const Vector edelta(-sin(delta)*cos(alpha), -sin(delta)*sin(alpha), cos(delta));

            Matrix I_inv_prime;
            I_inv_prime << -2.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 0.0,
                           -1.0/Ixy*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 2.0/Ixy*sin(delta)/cos(delta)/cos(delta)/cos(delta), 0.0,
                            0.0, 0.0, 0.0;

            for (unsigned j = 0; j < N; j++) {
                if (j == i) {continue;}
                const Vector& y = bodies.Getx(j);
                const double& My = bodies.GetM(j);

                r_vec = y - x;
                r = r_vec.norm();
                e_r = r_vec /r;
                r3_inv = 1/(r*r*r);
                ewer = ew.dot(e_r);

                torque(0) = 0.0;
                torque(1) = 3*PHYS_G*Mx*My*r3_inv*R*R*J2*ewer*(ealpha.dot(e_r))*cos(delta);
                torque(2) = -1.0/2*(I_inv_prime *L).dot(L)
                            + 3*PHYS_G*Mx*My*r3_inv*R*R*J2*ewer*(edelta.dot(e_r));

                Torque.row(i) += torque;
            }
        }

        //std::chrono::_V2::system_clock::time_point stop = std::chrono::high_resolution_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;

        return Torque;
    };


    // initialize solver
    Forest_Ruth<N> solver(x_dot, p_dot, orient_dot, L_dot, initial_bodies, t0);

    const double &t = solver.GetCurrentTime();
    const Body &sun = solver.GetCurrentBody(0);
    const Body &earth = solver.GetCurrentBody(1);
    const Body &moon = solver.GetCurrentBody(2);

    fs::path filename;
    bool solar_eclipse = false;
    bool solar_eclipse_active = false;
    double t_max = 0.0;
    double dist;
    double min_dist;
    std::map<std::string, std::vector<double>> general_data;
    std::array<std::vector<double>, 6> shadow_intersects;
    std::array<std::vector<u_int64_t>, 3> offsets;


    std::cout << "Starting integration with dt = " << dt << " seconds for a total time of " << T << " seconds (" << T/60/60/24/365.25 << " years) (" << K << " steps)." << "\n" << std::endl;

    // start integration loop
    for (size_t i = 0; i < K; i++) {
        solver.MakeStep(dt);

        solar_eclipse = eclipsed(sun, moon, earth);

        if (solar_eclipse && !solar_eclipse_active) {
            std::cout << "detected solar eclipse start at time " << j2000_to_utc_datetime(t) << "." << std::endl;

            filename = "solar_eclipse_of_" + j2000_to_utc_date(t_max);
            for (unsigned i : {0, 1, 2}) { // start offsets with 0
                offsets[i].push_back(0);
            }

            min_dist = INFINITY;
        }
        if (!solar_eclipse && solar_eclipse_active) {
            std::cout << "detected solar eclipse maximum at time " << j2000_to_utc_datetime(t_max) << "." << std::endl;
            std::cout << "detected solar eclipse end at time " << j2000_to_utc_datetime(t) << "." << "\n" << std::endl;

            write_eclipse_to_NetCDF(folder / filename,
                                    general_data,
                                    shadow_intersects,
                                    offsets);

            general_data = {};
            shadow_intersects = {};
            offsets = {};
        }
        solar_eclipse_active = solar_eclipse;

        if (solar_eclipse) {
            // compute distance of sun-moon axis to earths center to check for maximum eclipse time
            dist = compute_shadow_axis_distance_to_earth_center(earth, moon, sun);
            if (dist < min_dist) {
                min_dist = dist;
                t_max = t;
            }

            // save time and suns RA and Dec
            std::array<double, 3> sun_spherical= compute_spherical_seen_from_earth(sun, earth);
            std::array<double, 3> moon_spherical= compute_spherical_seen_from_earth(moon, earth);
            general_data["time"].push_back(t);
            general_data["sun_r"].push_back(sun_spherical[0]);
            general_data["sun_RA"].push_back(sun_spherical[1]);
            general_data["sun_Dec"].push_back(sun_spherical[2]);
            general_data["moon_r"].push_back(moon_spherical[0]);
            general_data["moon_RA"].push_back(moon_spherical[1]);
            general_data["moon_Dec"].push_back(moon_spherical[2]);

            // compute shadow intersections with earth and save
            std::array<std::vector<double>, 6> new_shadow_intersects = compute_shadow_earth_intersection(earth, moon, sun);
            for (unsigned i : {0, 1, 2, 3, 4, 5}) {
                shadow_intersects[i].insert(shadow_intersects[i].end(), new_shadow_intersects[i].begin(), new_shadow_intersects[i].end());
            }
            for (unsigned i : {0, 1, 2}) {
                offsets[i].push_back(shadow_intersects[2*i].size());
            }
        }




    }

    std::cout << "Integration complete." << std::endl;
}



