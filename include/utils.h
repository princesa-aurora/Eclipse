#include <array>
#include <stdexcept>
#include <cmath>
#include <string>

#include <erfa.h>
#include <Eigen/Dense>

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



double intersection_area(double R1, double R2, double d) {
    // return the area of intersection of two circles with radii R1 and R2, separated by a distance d
    if (d >= R1 + R2) {return 0.0;}

    double cos_phi = (R1*R1 - R2*R2 + d*d)/(2*d*R1);
    double cos_theta = (R2*R2 - R1*R1 + d*d)/(2*d*R2);

    double phi = acos(cos_phi);
    double theta = acos(cos_theta);

    double A1 = R1*R1*(phi - sin(phi)*cos_phi);
    double A2 = R2*R2*(theta - sin(theta)*cos_theta);

    return A1+A2;
}



double eclipse_ratio(const Vector &x0, double R0, const Vector &x1, double R1, const Vector &x2, double R2) {
    // check if body 1 ecplipses between bodies 0 and 2
    // return the occultation ratio of body 0 as seen from infinity
    // bad function, need to account for persepctive

    Vector e_parallel = (x2-x0)/(x2-x0).norm();

    double d01 = (x1-x0).dot(e_parallel);

    if (d01 < 0 || d01 > (x2-x0).norm()) {
        return false;
    }

    double D = (x1-x0 - d01*e_parallel).norm();
    double r_lightcone_1 = R0 + (R2-R0)*d01/(x2-x0).norm();

    double A_intersect = intersection_area(R1, r_lightcone_1, D);

    return A_intersect / (M_PI*R0*R0);
}



bool eclipsed(const Vector &x0, double R0, const Vector &x1, double R1, const Vector &x2, double R2) {
    // check if body 1 ecplipses between bodies 0 and 2

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


std::string jdtt_to_utc_date(double jdtt) {
    // convert Modified Julian Date (days) in TTT to a UTC date string

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
    return std::to_string(day) + "." + std::to_string(month) + "." + std::to_string(year) + " " +
           std::to_string(hour) + ":" + std::to_string(minute) + ":" + std::to_string(sec) + " UTC";
}

std::string j2000_to_utc_date(double j2000_seconds) {
    // convert seconds since J2000 to a UTC date string

    double jdtt = j2000_seconds/86400.0 + 2451545.0; // convert to Julian Days

    return jdtt_to_utc_date(jdtt);
}

double sun_zenith(double RA, double Dec, Body earth, Body sun) {
    // at RA, Dec on earths surface what is the zenith of the sun?
    Vector pos_on_earth(cos(Dec)*cos(RA), cos(Dec)*sin(RA), sin(Dec));
    Matrix RotMat = earth.GetRotMat();
    pos_on_earth = RotMat *pos_on_earth; // rotate from terrestrial coordinates to celestial ones

    Vector earth_to_sun = sun.Getx() - earth.Getx();
    double cos_zen = pos_on_earth.dot(earth_to_sun/earth_to_sun.norm());

    return M_PI/2 - acos(cos_zen);
}






#endif //UTILS_H_INCLUDED

