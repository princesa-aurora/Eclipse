#include <iostream>
#include <functional>
#include <vector>
#include <chrono>
#include <filesystem>
#include <format>
#include <cstdio>

#include <include/utils.h>
#include <include/initial_conditions.h>
#include <include/solver.h>

namespace fs = std::filesystem;


const size_t N = 10; // number of bodies in the system

const double t0 = 0.0; // seconds since J2000.0 epoch
const BodyArray<N> initial_bodies(heap_array{{Sun, Earth, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune}});

// start and end dates of the eclipse analysis season
const std::string start_date = "09.08.2026";
const std::string end_date = "15.08.2026";
const double T_start = utc_date_and_time_to_j2000(start_date, "00-00-00");
const double T_end = utc_date_and_time_to_j2000(end_date, "00-00-00");

// step sizes: dt1 default, dt2 during eclipse analysis
const double dt1 = 10.0;
const double dt2 = 1.0;

// number of points in Fibonacci sphere for eclipse analysis
const size_t grid_size = 5e6;

// refresh intervals for progress indicators
size_t progress_refresh_interval_catchup = 3000;
size_t progress_refresh_interval_analysis = 10;



int main() {
    // name and create output folder using a timestamp
    std::string timestamp = std::format("{:%d.%m.%Y_%H-%M-%S}",
                                        std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::system_clock::now()));
    fs::path folder = "/data/eclipse_data" / fs::path(timestamp);
    fs::create_directory(folder);


    std::function<VectorArray<N>(const BodyArray<N>&)> x_dot = [&](const BodyArray<N> &bodies)
    {
        // position time derivative
        VectorArray<N> v;
        for (unsigned i = 0; i < N; i++) {
            v.row(i) = bodies[i].Getv();
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
            const Vector& x = bodies[i].Getx();
            const double& Mx = bodies[i].GetM();
            const double& Rx = bodies[i].GetR();
            const double& J2x = bodies[i].GetJ2();
            const Vector ewx = bodies[i].GetAxis();

            for (unsigned j = 0; j < i; j++) {
                const Vector& y = bodies[j].Getx();
                const double& My = bodies[j].GetM();
                const double& Ry = bodies[j].GetR();
                const double& J2y = bodies[j].GetJ2();
                const Vector ewy = bodies[j].GetAxis();

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
            w.row(i) = bodies[i].Getw();
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
            const Vector& x = bodies[i].Getx();
            const double& alpha = bodies[i].Getorient()(1);
            const double& delta = bodies[i].Getorient()(2);
            const Vector& L = bodies[i].GetL();
            const double& Mx = bodies[i].GetM();
            const double& Ixy = bodies[i].GetIxy();
            const double& Iz = bodies[i].GetIz();
            const double& R = bodies[i].GetR();
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

    // create sampling grid for eclipse occultation data (Fibonacci sphere)
    heap_array<double, grid_size> lon_grid;
    heap_array<double, grid_size> lat_grid;
    for (size_t i = 0; i < grid_size; i++) {
        lon_grid[i] = fmod((3-sqrt(5))*M_PI *i, 2*M_PI);
        lat_grid[i] = asin(1.0 - double(2*i+1)/grid_size);
    }

    size_t iteration_counter = 0;
    bool solar_eclipse = false;
    bool analysis_active = false;
    double t_max;
    double dist;
    double min_dist;
    fs::path file_path;
    fs::path temp_path = folder / "eclipse_temp.nc";
    Eclipse_NetCDF<grid_size, 7> netcdf_file;
    std::array<std::string, 7> general_data_keys{"time", "r_sun", "lon_sun", "lat_sun", "r_moon", "lon_moon", "lat_moon"};
    std::array<double, 7> general_buffer;
    heap_array<double, grid_size> occult_buffer;
    heap_array<u_int8_t, grid_size> classif_buffer;

    std::cout << "Catching up to eclipse anaylsis window which starts at " << start_date << " 00-00-00 UTC." << std::endl;
    // do catching-up integration
    while (t <= T_start - dt1) {
        // keep going until just before the analysis window opens

        iteration_counter++;
        solver.MakeStep(dt1);

        // progress indicator
        if ((iteration_counter % progress_refresh_interval_catchup) == 0) {std::cout << "\r\033[2K" << "Currently at: " << j2000_to_utc_date(t) << std::flush;}
    }
    std::cout << "\r\033[2K" << std::endl;
    // if at the beginning of the analysis window an eclipse is ongoing skip ahead till it has ended
    // since that eclipse would be only recorded incompletely so we may just throw it away entirely
    if (eclipsed(earth, moon, sun)) {
        std::cout << "Note: at the beginning of the eclipse analysis window there is an ongoing eclipse."
                  << "As it can't be analyzed completely it is skipped entirely."
                  << "Did you maybe set the window slighly too late?" << std::endl;
        while (eclipsed(earth, moon, sun)) {
            solver.MakeStep(dt1);
        }
    }

    std::cout << "Reached the start of the eclipse analysis window.\n"
              << "Now checking for and analyzing eclipses until the end of the analysis window at " << end_date << " 00-00-00 UTC.\n" << std::endl;
    // do eclipse search and anaylsis
    while ((t < T_end) || analysis_active) {
        // keep going until the analysis window closes,
        // once the window is closed finish the last analysis before stopping

        iteration_counter++;
        solver.MakeStep(analysis_active ? dt2 : dt1);

        solar_eclipse = eclipsed(earth, moon, sun); // check for an eclipse

        if (solar_eclipse && !analysis_active) {
            // a solar eclipse is detected and we're currently not in analysis mode
            analysis_active = true; // start analysis mode

            std::cout << "\r\033[2K" << "Solar eclipse analysis triggered at " << j2000_to_utc_datetime(t) << ".\n"
                      << "Refining last step to narrow down start time." << std::endl;
            while(eclipsed(earth, moon, sun)) {
                solver.MakeStep(-dt2);
            }
            std::cout << "Detected solar eclipse start at time " << j2000_to_utc_datetime(t+dt2) << "." << std::endl;

            min_dist = INFINITY; // start min_dist at infinity
            netcdf_file.create_new_file(temp_path, general_data_keys, lon_grid, lat_grid); // create new file for the analysis
        }

        if (analysis_active) { // collection of analyses to be performed in analysis mode
            // compute distance of sun-moon axis to earths center to check for maximum eclipse time
            dist = compute_shadow_axis_distance_to_earth_center(earth, moon, sun);
            if (dist < min_dist) {
                min_dist = dist;
                t_max = t;
            }

            // save time and sun and moon positions
            std::array<double, 3> sun_spherical= compute_spherical_seen_from_earth(sun, earth);
            std::array<double, 3> moon_spherical= compute_spherical_seen_from_earth(moon, earth);
            general_buffer = {t,
                            sun_spherical[0], sun_spherical[1], sun_spherical[2],
                            moon_spherical[0], moon_spherical[1], moon_spherical[2]};

            // compute local occultations and save
            compute_local_occultations(earth, moon, sun, lon_grid, lat_grid,
                                                occult_buffer, classif_buffer);

            // write data to the file
            netcdf_file.write_step(general_buffer, occult_buffer, classif_buffer);
        }

        if (!solar_eclipse && analysis_active) {
            // no solar eclipse is detected and we're currently in analysis mode
            analysis_active = false; // else end analysis mode

            std::cout << "\r\033[2K" << "Detected solar eclipse maximum at time " << j2000_to_utc_datetime(t_max) << ".\n"
                      << "Detected solar eclipse end at time " << j2000_to_utc_datetime(t) << "." << "\n" << std::endl;

            netcdf_file.close_file(); // close output file
            file_path = folder / ("solar_eclipse_of_" + j2000_to_utc_date(t_max) + ".nc");
            // move output file from temp_path to file_path
            fs::rename(temp_path.c_str(), file_path.c_str());
        }

        // progress indicator
        if ((iteration_counter % progress_refresh_interval_analysis) == 0) {std::cout << "\r\033[2K" << "Currently at: " << j2000_to_utc_datetime(t) << std::flush;}
    }
    std::cout << "\r\033[2K" << std::flush;

    std::cout << "Reached the end of the eclipse analysis window.\n"
              << "Done, did a total of " << iteration_counter << " iterations." << std::endl;
}



