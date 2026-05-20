#include <array>
#include <stdexcept>
#include <cmath>
#include <string>
#include <vector>
#include <map>

#include <erfa.h>
#include <Eigen/Dense>
#include <netcdf.h>

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED


double PHYS_G = 6.6743015E+04; // km3/(1e24 kg)/s^2, Newtons constant of gravity


using Vector = Eigen::Vector3d;
using Matrix = Eigen::Matrix3d;

template<unsigned N>
using ScalarArray = Eigen::Vector<double, N>;
template<unsigned N>
using VectorArray = Eigen::Matrix<double, N, 3>;


class Body {

public:

    Body(
        double M,
        double a,
        double b,
        double i_f,
    
        Vector x0,
        Vector p0,
        Vector orient0,
        Vector L0) :
    M_(M),
    a_(a),
    b_(b),
    R_(pow(a*a*b, 1.0/3.0)),
    i_f_(i_f),
    Iz_(i_f*M*a*a),
    Ixy_(i_f/2.0*M*(a*a+b*b)),
    J2_(i_f/2.0*(a*a - b*b)/(R_*R_)),

    x_(x0),
    p_(p0),
    orient_(orient0),
    L_(L0)
    {}

    double GetM() const {
        return M_;
    }

    double Geta() const {
        return a_;
    }

    double Getb() const {
        return b_;
    }

    double GetR() const {
        return R_;
    }

    double Geti_f() const {
        return i_f_;
    }

    double GetIz() const {
        return Iz_;
    }

    double GetIxy() const {
        return Ixy_;
    }

    double GetJ2() const {
        return J2_;
    }


    const Vector& Getx() const {
        return x_;
    }

    const Vector& Getp() const {
        return p_;
    }

    const Vector& Getorient() const {
        return orient_;
    }

    const Vector& GetL() const {
        return L_;
    }


    Vector Getv() const {
        return p_ /M_;
    }

    Vector Getw() const {
        const double& delta = orient_(2);

        Matrix I_inv;
        I_inv << 1.0/Iz_ + 1.0/Ixy_*cos(delta)*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 1.0/Ixy_/(1.0+sin(delta)), 0.0,
                 1.0/Ixy_/(1.0+sin(delta)), 1.0/Ixy_/cos(delta)/cos(delta), 0.0,
                 0.0, 0.0, 1.0/Ixy_;

        return I_inv *L_;
    }


    Vector GetAxis() const {
        const double& alpha = orient_(1);
        const double& delta = orient_(2);

        return Vector(cos(delta)*cos(alpha), cos(delta)*sin(alpha), sin(delta));
    }

    double GetT_rot() const {
        double w = Getw()(0);

        return 2*M_PI /w;
    }

    Matrix GetRotMat() const {
        // get the rotation matrix that transforms from terrestrial coordinates to celestial ones
        const double& phi = orient_(0);
        const double& alpha = orient_(1);
        const double& delta = orient_(2);

        Matrix RotMat_phi;
        RotMat_phi << cos(phi), -sin(phi), 0.0,
                      sin(phi), cos(phi), 0.0,
                      0.0, 0.0, 1.0;
        Matrix RotMat_axis;
        RotMat_axis << 1 - cos(alpha)*cos(alpha)*(1.0-sin(delta)), -cos(alpha)*sin(alpha)*(1.0-sin(delta)), cos(delta)*cos(alpha),
                       -cos(alpha)*sin(alpha)*(1.0-sin(delta)), 1 - sin(alpha)*sin(alpha)*(1.0-sin(delta)), cos(delta)*sin(alpha),
                       -cos(delta)*cos(alpha), -cos(delta)*sin(alpha), sin(delta);

        return RotMat_axis*RotMat_phi;
    }


    void Setx(const Vector &x) {
        x_ = x;
    }

    void Incrementx(const Vector &incr_x) {
        x_ += incr_x;
    }

    void Setp(const Vector &p) {
        p_ = p;
    }

    void Incrementp(const Vector &incr_p) {
        p_ += incr_p;
    }

    void Setorient(const Vector &orient) {
        orient_ = orient;
        enforce_orient_range_();
    }

    void Incrementorient(const Vector &incr_orient) {
        orient_ += incr_orient;
        enforce_orient_range_();
    }

    void SetL(const Vector &L) {
        L_ = L;
    }

    void IncrementL(const Vector &incr_L) {
        L_ += incr_L;
    }


private:
    // properties of the body
    const double M_; // mass
    const double a_; // equatorial radius
    const double b_; // polar radius
    const double R_; // volumetric mean radius
    const double i_f_; // moment of inertia factor
    const double Iz_; // moment of inertia along polar axis
    const double Ixy_; // moment of inertia along equatorial axes
    const double J2_; // second zonal coefficient

    // dynamic variables of the body
    Vector x_; // position
    Vector p_; // translational momentum
    Vector orient_; // orientation: (0): rotation angle, (1): pole RA, (2): pole Dec
    Vector L_; // angular momentum: components corresponding to orientation vector

    void enforce_orient_range_() {
        double& alpha = orient_(1);
        double& delta = orient_(2);

        delta = fmod(delta, 2*M_PI);
        if (delta > 3*M_PI/2) { // 3pi/2 < delta < 2pi
            delta = delta - 2*M_PI;
        }
        else if (delta > M_PI/2) { // pi/2 < delta < 3pi/2
            delta = M_PI - delta;
            alpha = alpha + M_PI;
        }
        // else: 0 < delta < pi/2 : continue

        alpha = fmod(alpha, 2*M_PI);
    }
};


template<unsigned N>
class BodyArray {

public:
    BodyArray(std::array<Body ,N> body_arr) :
        body_arr_(body_arr)
    {}

    Body& operator[](unsigned idx) {
        return body_arr_[idx];
    }

    const Body& operator[](unsigned idx) const {
        return body_arr_[idx];
    }


    double GetM(unsigned idx) const {

        return body_arr_[idx].GetM();
    }

    double Geta(unsigned idx) const {

        return body_arr_[idx].Geta();
    }

    double Getb(unsigned idx) const {

        return body_arr_[idx].Getb();
    }

    double GetR(unsigned idx) const {

        return body_arr_[idx].GetR();
    }

    double Geti_f(unsigned idx) const {

        return body_arr_[idx].Geti_f();
    }

    double GetIz(unsigned idx) const {

        return body_arr_[idx].GetIz();
    }

    double GetIxy(unsigned idx) const {

        return body_arr_[idx].GetIxy();
    }

    double GetJ2(unsigned idx) const {

        return body_arr_[idx].GetJ2();
    }


    const Vector& Getx(unsigned idx) const {

        return body_arr_[idx].Getx();
    }

    const Vector& Getp(unsigned idx) const {

        return body_arr_[idx].Getp();
    }

    const Vector& Getorient(unsigned idx) const {

        return body_arr_[idx].Getorient();
    }

    const Vector& GetL(unsigned idx) const {

        return body_arr_[idx].GetL();
    }


    Vector Getv(unsigned idx) const {

        return body_arr_[idx].Getv();
    }

    Vector Getw(unsigned idx) const {

        return body_arr_[idx].Getw();
    }


    Vector GetAxis(unsigned idx) const {

        return body_arr_[idx].GetAxis();
    }

    double GetT_rot(unsigned idx) const {

        return body_arr_[idx].GetT_rot();
    }

    Matrix GetRotMat(unsigned idx) const {

        return body_arr_[idx].GetRotMat();
    }


    void Setx(const VectorArray<N> &x) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setx(x.row(i));
        }
    }

    void Incrementx(const VectorArray<N> &incr_x) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementx(incr_x.row(i));
        }
    }

    void Setp(const VectorArray<N> &p) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setp(p.row(i));
        }
    }

    void Incrementp(const VectorArray<N> &incr_p) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementp(incr_p.row(i));
        }
    }

    void Setorient(const VectorArray<N> &orient) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setorient(orient.row(i));
        }
    }

    void Incrementorient(const VectorArray<N> &incr_orient) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementorient(incr_orient.row(i));
        }
    }

    void SetL(const VectorArray<N> &L) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].SetL(L.row(i));
        }
    }

    void IncrementL(const VectorArray<N> &incr_L) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].IncrementL(incr_L.row(i));
        }
    }


private:
    std::array<Body, N> body_arr_;
};



bool eclipsed(const Body &body0, const Body &body1, const Body &body2) {
    // check if body 1 ecplipses between bodies 0 and 2
    const double &R0 = body0.GetR();
    const double &R1 = body1.GetR();
    const double &R2 = body2.GetR();
    const Vector &x0 = body0.Getx();
    const Vector &x1 = body1.Getx();
    const Vector &x2 = body2.Getx();

    Vector e_parallel = (x2-x0)/(x2-x0).norm();

    double d01 = (x1-x0).dot(e_parallel);

    if (d01 < 0 || d01 > (x2-x0).norm()) {
        return false;
    }

    double D = (x1-x0 - d01*e_parallel).norm();
    double r_lightcone_1 = R0 + (R2-R0)*d01/(x2-x0).norm();

    if (D < r_lightcone_1 + R1) {
        return true;
    }
    else {
        return false;
    }
}


std::pair<std::string, std::string> jdtt_to_utc_date_and_time(double jdtt) {
    // convert Modified Julian Date (days) in TT to a UTC date and a time string

    double tt_int = floor(jdtt); // split into integer
    double tt_frac = jdtt - tt_int; // and fractional parts

    // convert TT to TAI (International Atomic Time)
    double tai_int, tai_frac;
    eraTttai(tt_int, tt_frac, &tai_int, &tai_frac);

    // convert TAI to UTC (Universal Coordinated Time)
    double utc_int, utc_frac;
    eraTaiutc(tai_int, tai_frac, &utc_int, &utc_frac);

    // break down UTC into year, month, day and a fraction of a day
    int year, month, day;
    double day_fraction;
    eraJd2cal(utc_int, utc_frac, &year, &month, &day, &day_fraction);

    // convert the fraction of a day into hours, minutes and seconds
    int hour = day_fraction * 24;
    int minute = day_fraction * 24*60 - hour*60;
    int sec = day_fraction * 24*60*60 - hour*60*60 - minute*60;

    // format into date string and return
    std::string date = std::to_string(day) + "." + std::to_string(month) + "." + std::to_string(year);
    std::string time = std::to_string(hour) + ":" + std::to_string(minute) + ":" + std::to_string(sec);
    
    return std::make_pair(date, time);
}

std::string j2000_to_utc_datetime(double j2000_seconds) {
    // convert seconds since J2000 to a UTC datetime string

    double jdtt = j2000_seconds/86400.0 + 2451545.0; // convert to Julian Days
    std::pair<std::string, std::string> date_and_time = jdtt_to_utc_date_and_time(jdtt);

    return date_and_time.first + " " + date_and_time.second + " UTC";
}

std::string j2000_to_utc_date(double j2000_seconds) {
    // convert seconds since J2000 to a UTC date string

    double jdtt = j2000_seconds/86400.0 + 2451545.0; // convert to Julian Days
    std::pair<std::string, std::string> date_and_time = jdtt_to_utc_date_and_time(jdtt);

    return date_and_time.first;
}



double compute_sun_zenith(double RA, double Dec, const Body &earth, const Body &sun) {
    // at RA, Dec on earths surface what is the zenith of the sun?

    Vector pos_on_earth(cos(Dec)*cos(RA), cos(Dec)*sin(RA), sin(Dec));
    Matrix RotMat_earth = earth.GetRotMat();
    pos_on_earth = RotMat_earth *pos_on_earth; // rotate from terrestrial frame to ICRF

    Vector earth_to_sun = sun.Getx() - earth.Getx();
    double cos_zen = pos_on_earth.dot(earth_to_sun/earth_to_sun.norm());

    return M_PI/2 - acos(cos_zen);
}

std::array<double, 3> compute_spherical_seen_from_earth(const Body &body, const Body &earth) {
    // what is the bodies distance, RA and Dec as seen from the terrestrial frame?

    Vector r = body.Getx() - earth.Getx();
    Matrix RotMat_inv_earth = earth.GetRotMat().transpose();
    r = RotMat_inv_earth *r; // rotate from ICRF to terrestrial frame

    double dist = r.norm();
    double RA = atan2(r(1), r(0));
    double Dec = atan(r(2)/sqrt(r(0)*r(0) + r(1)*r(1)));

    return {dist, RA, Dec};
}

double compute_shadow_axis_distance_to_earth_center(const Body &earth, const Body &moon, const Body &sun) {
    // the shadow's central axis is given be the line connecting sun and moon
    // return the shortest distance of this line to the earths center

    Vector moon_to_earth = earth.Getx() - moon.Getx();
    Vector sun_to_moon = moon.Getx() - sun.Getx();
    sun_to_moon /= sun_to_moon.norm();

    double dist2 = moon_to_earth.dot(moon_to_earth) - (moon_to_earth.dot(sun_to_moon))*(moon_to_earth.dot(sun_to_moon));

    return sqrt(dist2);
}



std::array<std::vector<double>, 6> compute_shadow_earth_intersection(const Body &earth, const Body &moon,  const Body &sun) {
    // compute the intersection curves between the umbra, antumbra and penumbra and the earths surface,
    // save RA and DEC of the curve points in the terrestrial frame in netcdf
    const double& R_earth = earth.GetR();
    const double& R_moon = moon.GetR();
    const double& R_sun = sun.GetR();

    Vector x_earth = earth.Getx();
    Vector x_moon = moon.Getx();
    const Vector& x_earth_ICRF = earth.Getx();
    const Vector& x_sun_ICRF = sun.Getx();

    const Matrix RotMat_inv_earth = earth.GetRotMat().transpose();

    // offset origin to be at the sun
    x_earth -= x_sun_ICRF;
    x_moon -= x_sun_ICRF;
    // x_sun = 0

    double r_moon = x_moon.norm();
    double theta_moon = asin(x_moon(2)/r_moon);
    double phi_moon = atan2(x_moon(1), x_moon(0));

    Matrix R_z; // rotate frame to align moon_to_sun with the z-axis
    R_z << 1 - cos(phi_moon)*cos(phi_moon)*(1.0-sin(theta_moon)), -cos(phi_moon)*sin(phi_moon)*(1.0-sin(theta_moon)), -cos(theta_moon)*cos(phi_moon),
                       -cos(phi_moon)*sin(phi_moon)*(1.0-sin(theta_moon)), 1 - sin(phi_moon)*sin(phi_moon)*(1.0-sin(theta_moon)), -cos(theta_moon)*sin(phi_moon),
                       cos(theta_moon)*cos(phi_moon), cos(theta_moon)*sin(phi_moon), sin(theta_moon);
    x_earth = R_z*x_earth;
    // x_moon = r_moon * e_z
    // x_sun = 0

    // compute cylindrical form of x_earth
    double rho_earth = sqrt(x_earth(0)*x_earth(0) + x_earth(1)*x_earth(1));
    double beta_earth = atan2(x_earth(1), x_earth(0));
    double z_earth = x_earth(2);

    // calculate umbral shadow parameters (shadow is a cone: r(t, beta) = t*(A*cos(beta), A*sin(beta), 1))+(0, 0, z0)
    // antumbra is just the continuation of the umbra beyond the vertex, i.e. no separate parameters are needed
    double z0_umbral = R_sun*r_moon/(R_sun - R_moon);
    double A_umbral = (R_sun-R_moon)/sqrt(r_moon*r_moon - (R_sun-R_moon)*(R_sun-R_moon));
    double t_moon_umbral = (r_moon-z0_umbral)/(A_umbral*A_umbral+1.0); // t at umbral shadow intersection with the moon

    // calculate penumbral shadow parameters
    double z0_penumbral = R_sun*r_moon/(R_sun+R_moon);
    double A_penumbral = (R_sun+R_moon)/sqrt(r_moon*r_moon - (R_sun+R_moon)*(R_sun+R_moon));
    double t_moon_penumbral = (r_moon-z0_penumbral)/(A_penumbral*A_penumbral+1.0); // t at penumbral shadow intersection with the moon

    // sample from intersection of shadow cones with earth: r(t(beta), beta)
    // thereby exclude points located on the side of earth that is facing away from the sun (it's dark there anyways)
    // umbral and antumbral intersections are split as antumbral: (t(beta) > 0) and umbral: (t_umbral_moon < t(beta) < 0)
    std::vector<double> umbral_RA;
    std::vector<double> umbral_Dec;
    std::vector<double> antumbral_RA;
    std::vector<double> antumbral_Dec;
    std::vector<double> penumbral_RA;
    std::vector<double> penumbral_Dec;
    double beta;
    double a_mn, b_mn, c_mn, D_mn;
    Vector x;
    double RA, Dec;

    unsigned K = 2000;
    for (unsigned k = 0; k < K; k++) {
        beta = 2*M_PI*k/K;

        // start with umbral/antumbral shadow
        a_mn = A_umbral*A_umbral + 1.0;
        b_mn = 2*(z0_umbral-z_earth - rho_earth*A_umbral*cos(beta-beta_earth));
        c_mn = (z0_umbral-z_earth)*(z0_umbral-z_earth) + rho_earth*rho_earth - R_earth*R_earth;
        D_mn = b_mn*b_mn - 4*a_mn*c_mn;

        if (D_mn >= 0.0) { // D_mn is positive => intersection exists

            // consider the positive and negative branches of the midnight formula
            for (double t : {(-b_mn + sqrt(D_mn))/2/a_mn, (-b_mn - sqrt(D_mn))/2/a_mn}) {
                x = Vector(t*A_umbral*cos(beta), t*A_umbral*sin(beta), t+z0_umbral);

                if ((x-x_earth).dot(x_earth) < 0.0) { // facing toward the sun
                    // undo all the previous transformations to get x in ICRF
                    x = R_z.transpose()*x;
                    x += x_sun_ICRF;
                    // apply the inverse of earths rotation matrix to transform into the terrestrial frame
                    x = RotMat_inv_earth*(x-x_earth_ICRF);

                    // now extract RA and Dec
                    RA = atan2(x(1), x(0));
                    Dec = atan(x(2)/sqrt(x(0)*x(0)+x(1)*x(1)));

                    // save the results while distinguishing between antumbral (t > 0) and umbral (t_moon < t < 0) intersections,
                    // if t < t_moon_umbral no eclipe happens, so don't save
                    if (t > 0.0) {
                        antumbral_RA.push_back(RA);
                        antumbral_Dec.push_back(Dec);
                    }
                    else if (t > t_moon_umbral) {
                        umbral_RA.push_back(RA);
                        umbral_Dec.push_back(Dec);
                    }
                }
            }
        }

        // now do the penumbral shadow
        a_mn = A_penumbral*A_penumbral + 1.0;
        b_mn = 2*(z0_penumbral-z_earth - rho_earth*A_penumbral*cos(beta-beta_earth));
        c_mn = (z0_penumbral-z_earth)*(z0_penumbral-z_earth) + rho_earth*rho_earth - R_earth*R_earth;
        D_mn = b_mn*b_mn - 4*a_mn*c_mn;

        if (D_mn >= 0.0) { // D_mn is positive => intersection exists

            // consider the positive and negative branches of the midnight formula
            for (double t : {(-b_mn + sqrt(D_mn))/2/a_mn, (-b_mn - sqrt(D_mn))/2/a_mn}) {
                x = Vector(t*A_penumbral*cos(beta), t*A_penumbral*sin(beta), t+z0_penumbral);

                if ((x-x_earth).dot(x_earth) < 0.0) { // facing toward the sun
                    // undo all the previous transformations to get x in ICRF
                    x = R_z.transpose()*x;
                    x += x_sun_ICRF;
                    // apply the inverse of earths rotation matrix to transform into the terrestrial frame
                    x = RotMat_inv_earth*(x-x_earth_ICRF);

                    // now extract RA and Dec
                    RA = atan2(x(1), x(0));
                    Dec = atan(x(2)/sqrt(x(0)*x(0)+x(1)*x(1)));

                    // save the results
                    // if t < t_moon_penumbral no eclipe happens, so don't save
                    if (t > t_moon_penumbral) {
                        penumbral_RA.push_back(RA);
                        penumbral_Dec.push_back(Dec);
                    }
                }
            }
        }
    }

    return {umbral_RA, umbral_Dec, antumbral_RA, antumbral_Dec, penumbral_RA, penumbral_Dec};
}



void write_eclipse_to_NetCDF(std::string file_path,
                             const std::map<std::string, std::vector<double>> &general_data,
                             const std::array<std::vector<double>, 6> &points,
                             const std::array<std::vector<u_int64_t>, 3> &offsets)
{
    // create NetCDF file
    int ncid;
    bool status = nc_create(file_path.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid);
    if (status != NC_NOERR) {
        throw std::runtime_error("NetCDF: Failed to create file: " + file_path);
    }

    if (general_data.size() == 0) {
        throw std::invalid_argument("NetCDF: general_data is empty, but it should at least contain time information.");
    }

    size_t num_steps = general_data.begin()->second.size();
    for (const std::pair<std::string, std::vector<double>> &data : general_data) {
        if ((data.second.size() != num_steps)) {
            throw std::invalid_argument("NetCDF: " +
                                        data.first + "(" + std::to_string(data.second.size()) + ") must have the same size as " +
                                        general_data.begin()->first + "(" + std::to_string(num_steps) + ").");
        }
    }

    // define general dimension and variables
    int general_dim;
    std::map<std::string, int> general_ids_map;
    nc_def_dim(ncid, "general", num_steps, &general_dim);

    for (const std::pair<std::string, std::vector<double>> data : general_data) {
        nc_def_var(ncid, data.first.c_str(), NC_DOUBLE, 1, &general_dim, &general_ids_map[data.first]);
    }

    // create groups for each shadow type and define variables in them
    std::array<std::string, 3> grp_names{"umbral", "antumbral", "penumbral"};
    std::array<std::array<int, 6>, 3> grp_ids;

    for (unsigned i : {0, 1, 2}) {
        std::string grp_name = grp_names[i];
        int& grp_id = grp_ids[i][0];
        int& RA_id = grp_ids[i][1];
        int& Dec_id = grp_ids[i][2];
        int& offsets_id = grp_ids[i][3];
        int& points_dim = grp_ids[i][4];
        int& offsets_dim = grp_ids[i][5];

        const std::vector<double> &RA = points[2*i];
        const std::vector<double> &Dec = points [2*i+1];
        const std::vector<u_int64_t> &offset = offsets[i];
        if (RA.size() != Dec.size()) { // check that RA and Dec have same length
            throw std::invalid_argument("NetCDF: group " + grp_name + ": RA(" + std::to_string(RA.size()) + ") and Dec(" + std::to_string(Dec.size()) + ") must have the same size.");
        }
        if (offset.size() != num_steps +1) {
            throw std::invalid_argument("NetCDF: group " + grp_name + ": offset(" + std::to_string(offset.size()) + ") must be one longer than times(" + std::to_string(num_steps) + ").");
        }

        nc_def_grp(ncid, grp_name.c_str(), &grp_id);

        nc_def_dim(grp_id, "points", Dec.size(), &points_dim);
        nc_def_dim(grp_id, "offsets", offset.size(), &offsets_dim);

        nc_def_var(grp_id, "RA", NC_DOUBLE, 1, &points_dim, &RA_id);
        nc_def_var(grp_id, "Dec", NC_DOUBLE, 1, &points_dim, &Dec_id);
        nc_def_var(grp_id, "offsets", NC_UINT64, 1, &offsets_dim, &offsets_id);
    }

    // close file definition mode
    nc_enddef(ncid);

    // write the general data
    for (const std::pair<std::string, std::vector<double>> data : general_data) {
        nc_put_var_double(ncid, general_ids_map[data.first], data.second.data());
    }

    // write data in each group
    for (unsigned i : {0, 1, 2}) {
        std::string grp_name = grp_names[i];
        int& grp_id = grp_ids[i][0];
        int& RA_id = grp_ids[i][1];
        int& Dec_id = grp_ids[i][2];
        int& offsets_id = grp_ids[i][3];

        const std::vector<double> &RA = points[2*i];
        const std::vector<double> &Dec = points [2*i+1];
        const std::vector<u_int64_t> &offset = offsets[i];

        // write RA and Dec
        nc_put_var_double(grp_id, RA_id, RA.data());
        nc_put_var_double(grp_id, Dec_id, Dec.data());

        // write offsets
        nc_put_var_ulonglong(grp_id, offsets_id, reinterpret_cast<const unsigned long long*>(offset.data()));
    }

    // close file
    if (ncid >= 0) {nc_close(ncid);}
}






#endif //UTILS_H_INCLUDED

