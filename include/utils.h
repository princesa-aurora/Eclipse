#include <iostream>
#include <array>
#include <stdexcept>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <cstdint>
#include <sstream>
#include <algorithm>

#include <erfa.h>
#include <Eigen/Dense>
#include <netcdf.h>

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED



using Vector = Eigen::Vector3d;
using Matrix = Eigen::Matrix3d;

class Quaternion : public Eigen::Vector4d {
public:
    typedef Eigen::Vector4d Base;

    Quaternion()
        : Base() {}

    Quaternion(const Base& data)
        : Base(data) {}

    Quaternion(double scal, double x, double y, double z)
        : Base(scal, x, y, z) {}

    Quaternion(double scal, const Vector& vec)
        : Base(scal, vec.x(), vec.y(), vec.z()) {}


    static Quaternion Real(double scal) {
        return Quaternion(scal, Vector::Zero());
    }

    static Quaternion Pure(const Vector& vec) {
        return Quaternion(0.0, vec);
    }

    static Quaternion Constant(double val) {
        return Quaternion(val, val, val, val);
    }

    static Quaternion Zero() {
        return Quaternion::Constant(0.0);
    }

    static Quaternion Identity() {
        return Quaternion(1.0, 0.0, 0.0, 0.0);
    }

    static Quaternion i() {
        return Quaternion(0.0, 1.0, 0.0, 0.0);
    }

    static Quaternion j() {
        return Quaternion(0.0, 0.0, 1.0, 0.0);
    }

    static Quaternion k() {
        return Quaternion(0.0, 0.0, 0.0, 1.0);
    }


    template<typename OtherDerived>
    Quaternion(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}

    template<typename OtherDerived>
    Quaternion& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        this->Base::operator=(other);
        return *this;
    }


    double& scalar() {
        return (*this)(0); 
    }

    const double& scalar() const {
        return (*this)(0);
    }

    Eigen::Ref<Vector> vector() {
        return (*this).tail<3>();
    }

    Eigen::Ref<const Vector> vector() const {
        return (*this).tail<3>();
    }


    Quaternion operator+(const Quaternion& other) const {
        return Quaternion(this->Base::operator+(other));
    }

    Quaternion& operator+=(const Quaternion& other) {
        this->Base::operator+=(other);
        return *this;
    }

    Quaternion operator-(const Quaternion& other) const {
        return Quaternion(this->Base::operator-(other));
    }

    Quaternion& operator-=(const Quaternion& other) {
        this->Base::operator-=(other);
        return *this;
    }

    Quaternion operator*(double scalar) const {
        return Quaternion(this->Base::operator*(scalar));
    }

    Quaternion& operator*=(double scalar) {
        this->Base::operator*=(scalar);
        return *this;
    }

    Quaternion operator/(double scalar) const {
        return Quaternion(this->Base::operator/(scalar));
    }

    Quaternion& operator/=(double scalar) {
        this->Base::operator/=(scalar);
        return *this;
    }


    Quaternion conjugate() const {
        return Quaternion(this->scalar(), -this->vector());
    }

    Quaternion inverse() const {
        double squared_norm = this->squaredNorm();
        return this->conjugate() /squared_norm;
    }

    Quaternion operator*(const Quaternion& other) const {
        Quaternion result;

        result.scalar() = this->scalar()*other.scalar() - this->vector().dot(other.vector());
        result.vector() = this->scalar()*other.vector() + this->vector()*other.scalar() + this->vector().cross(other.vector());

        return result;
    }

    Quaternion& operator*=(const Quaternion& other) {
        *this = *this * other;
        return *this;
    }

    Quaternion operator*(const Vector& vec) const {
        Quaternion result;

        result.scalar() = -this->vector().dot(vec);
        result.vector() = this->scalar()*vec + this->vector().cross(vec);

        return result;
    }

    Quaternion& operator*=(const Vector& vec) {
        *this = *this * vec;
        return *this;
    }

    Vector rotate(const Vector& vec) const {
        return ((*this) * vec * this->conjugate()).vector();
    }

    Quaternion mult_i() const {
        Quaternion result;

        result.scalar() = -this->vector().x();
        result.vector().x() = this->scalar();
        result.vector().y() = this->vector().z();
        result.vector().z() = -this->vector().y();

        return result;
    }

    Quaternion mult_j() const {
        Quaternion result;

        result.scalar() = -this->vector().y();
        result.vector().x() = -this->vector().z();
        result.vector().y() = this->scalar();
        result.vector().z() = this->vector().x();

        return result;
    }

    Quaternion mult_k() const {
        Quaternion result;

        result.scalar() = -this->vector().z();
        result.vector().x() = this->vector().y();
        result.vector().y() = -this->vector().x();
        result.vector().z() = this->scalar();

        return result;
    }
};

Quaternion operator*(double scalar, const Quaternion& q) {
    return q *scalar;
}

Quaternion operator*(const Vector& vec, const Quaternion& q) {
    Quaternion result;

    result.scalar() = -vec.dot(q.vector());
    result.vector() = vec*q.scalar() + vec.cross(q.vector());

    return result;
}


template<unsigned N>
using ScalarArray = Eigen::Matrix<double, N, 1>;
template<unsigned N>
using VectorArray = Eigen::Matrix<double, N, 3>;
template<unsigned N>
using QuaternionArray = Eigen::Matrix<double, N, 4>;



std::pair<Quaternion, Quaternion> angles_to_quaternions(double phi, double alpha, double delta,
                                double phi_dot, double alpha_dot, double delta_dot) {
    // represent the orientation angles as a quaternion (and same for the time derivative)
    Vector axis(cos(delta)*cos(alpha), cos(delta)*sin(alpha), sin(delta));
    Vector axis_dot(-sin(delta)*delta_dot*cos(alpha) -cos(delta)*sin(alpha)*alpha_dot,
                    -sin(delta)*delta_dot*sin(alpha) +cos(delta)*cos(alpha)*alpha_dot,
                     cos(delta)*delta_dot);

    Quaternion q(cos(phi/2), sin(phi/2) *axis);

    Quaternion w(-sin(phi/2) *phi_dot/2, cos(phi/2) *phi_dot/2 *axis + sin(phi/2) *axis_dot);

    return {q, w};
}



template <typename T, size_t k>
class heap_array { // basically like a std::array, but allocates memory on the heap not stack for much larger possible sizes
public:
    // default constructor
    heap_array() : data_(k) {} // allocate memory for the desired k elements

    // value assignment constructor
    heap_array(const T (&init_array)[k]) : data_(init_array, init_array + k) {}

    // expose standard array methods
    T& operator[](size_t index) {
        return data_[index];
    }
    const T& operator[](size_t index) const {
        return data_[index];
    }
    T* data() {
        return data_.data();
    }
    const T* data() const {
        return data_.data();
    }
    size_t size() const {
        return k;
    }

private:
    std::vector<T> data_; // std::vector allocates memory on the heap
};



class Body {

public:

    Body(
        std::string name,

        double M,
        double a,
        double b,
        double i_f,
        double J2,

        Vector x0,
        Vector v0,
        Quaternion q0,
        Quaternion w0) :
    name_(name),
    M_(M),
    a_(a),
    b_(b),
    R_(pow(a*a*b, 1.0/3.0)),
    i_f_(i_f),
    Iz_(i_f*M*a*a),
    Ixy_(i_f/2.0*M*(a*a+b*b)),
    J2_(0.0), //J2),

    x_(x0),
    p_(Vector::Constant(NAN)),
    v_(v0),
    q_(q0),
    L_(Quaternion::Constant(NAN)),
    w_(w0)
    {}

    std::string Getname() const {
        return name_;
    }

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

    const Quaternion& Getq() const {
        return q_;
    }

    const Quaternion& GetL() const {
        return L_;
    }


    const Vector& Getv() const {
        return v_;
    }

    const Quaternion& Getw() const {
        return w_;
    }


    Vector GetPoleAxis() const {
        Vector e_z = Vector(0.0, 0.0, 1.0);
        return q_.rotate(e_z);
    }

    Vector Getw_phys() const {
        return 2*(w_*q_.conjugate()).vector();
    }

    double GetT_rot() const {
        double omega = 2*w_.norm();

        return 2*M_PI /omega;
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

    void Setv(const Vector &v) {
        v_ = v;
    }

    void Incrementv(const Vector &incr_v) {
        v_ += incr_v;
    }

    void Setq(const Quaternion &q) {
        q_ = q;
        //q_.normalize(); // enforce unit quaternion (actually screws things up)
    }

    void Incrementq(const Quaternion &incr_q) {
        q_ += incr_q;
        //q_.normalize(); // enforce unit quaternion (actually screws things up)
    }

    void SetL(const Quaternion &L) {
        L_ = L;
    }

    void IncrementL(const Quaternion &incr_L) {
        L_ += incr_L;
    }

    void Setw(const Quaternion &w) {
        w_ = w;
    }

    void Incrementw(const Quaternion &incr_w) {
        w_ += incr_w;
    }


private:
    // name of the body
    const std::string name_;

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
    Vector v_; // velocity
    Quaternion q_; // orientation quaternion
    Quaternion L_; // angular momentum (conjugate momentum to q)
    Quaternion w_; // angular velocity (time derivative of q)
};


template<unsigned N>
class BodyArray {

public:
    BodyArray(heap_array<Body ,N> body_arr) :
        body_arr_(body_arr)
    {}

    Body& operator[](unsigned idx) {
        return body_arr_[idx];
    }

    const Body& operator[](unsigned idx) const {
        return body_arr_[idx];
    }


    VectorArray<N> Getx() const {
        VectorArray<N> x_arr;
        for (unsigned i = 0; i < N; i++) {
            x_arr.row(i) = body_arr_[i].Getx();
        }
        return x_arr;
    }

    VectorArray<N> Getp() const {
        VectorArray<N> p_arr;
        for (unsigned i = 0; i < N; i++) {
            p_arr.row(i) = body_arr_[i].Getp();
        }
        return p_arr;
    }

    VectorArray<N> Getv() const {
        VectorArray<N> v_arr;
        for (unsigned i = 0; i < N; i++) {
            v_arr.row(i) = body_arr_[i].Getv();
        }
        return v_arr;
    }

    QuaternionArray<N> Getq() const {
        QuaternionArray<N> q_arr;
        for (unsigned i = 0; i < N; i++) {
            q_arr.row(i) = body_arr_[i].Getq();
        }
        return q_arr;
    }

    QuaternionArray<N> GetL() const {
        QuaternionArray<N> L_arr;
        for (unsigned i = 0; i < N; i++) {
            L_arr.row(i) = body_arr_[i].GetL();
        }
        return L_arr;
    }

    QuaternionArray<N> Getw() const {
        QuaternionArray<N> w_arr;
        for (unsigned i = 0; i < N; i++) {
            w_arr.row(i) = body_arr_[i].Getw();
        }
        return w_arr;
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

    void Setv(const VectorArray<N> &v) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setv(v.row(i));
        }
    }

    void Incrementv(const VectorArray<N> &incr_v) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementv(incr_v.row(i));
        }
    }

    void Setq(const QuaternionArray<N> &q) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setq(q.row(i));
        }
    }

    void Incrementq(const QuaternionArray<N> &incr_q) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementq(incr_q.row(i));
        }
    }

    void SetL(const QuaternionArray<N> &L) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].SetL(L.row(i));
        }
    }

    void IncrementL(const QuaternionArray<N> &incr_L) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].IncrementL(incr_L.row(i));
        }
    }

    void Setw(const QuaternionArray<N> &w) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setw(w.row(i));
        }
    }

    void Incrementw(const QuaternionArray<N> &incr_w) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementw(incr_w.row(i));
        }
    }


private:
    heap_array<Body, N> body_arr_;
};



bool eclipsed(const Body &earth, const Body &moon, const Body &sun) {
    // check if body 1 ecplipses between bodies 0 and 2
    const double& R_earth = earth.GetR();
    const double& R_moon = moon.GetR();
    const double& R_sun = sun.GetR();

    const Vector& x_earth_ICRF = earth.Getx();
    const Vector& x_moon_ICRF = moon.Getx();
    const Vector& x_sun_ICRF = sun.Getx();
    const Vector& v_moon = moon.Getv();
    const Vector& v_sun = sun.Getv();

    // shift ICRF origin to earth
    // x_earth = 0;
    Vector x_moon = x_moon_ICRF - x_earth_ICRF;
    Vector x_sun = x_sun_ICRF - x_earth_ICRF;

    // we need the retarded positions of the sun and moon (ie corrected for light travel time)
    double PHYS_c = 299792.458; // remove before merge!!!!
    double t_moon = -1.0/PHYS_c *x_moon.norm();
    double t_sun = -1.0/PHYS_c *x_sun.norm();
    x_moon += v_moon*t_moon;
    x_sun += v_sun*t_sun;

    double d_sun = x_sun.norm();
    Vector e_parallel = x_sun/d_sun;
    double d_moon_parallel = x_moon.dot(e_parallel);
    double d_moon_orthogonal = sqrt(x_moon.dot(x_moon) - d_moon_parallel*d_moon_parallel);

    double z0 = R_earth*d_sun/(R_earth - R_sun);
    double A = R_earth/sqrt(z0*z0 - R_earth*R_earth);

    double d_moon_orthogonal_boundary = sqrt(A*A+1)*R_moon + A*(d_moon_parallel - z0);

    return d_moon_orthogonal <= d_moon_orthogonal_boundary;
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
    std::string date = std::format("{:02}", day) + "." + std::format("{:02}", month) + "." + std::format("{:02}", year);
    std::string time = std::format("{:02}", hour) + "-" + std::format("{:02}", minute) + "-" + std::format("{:02}", sec);

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


double utc_date_and_time_to_jdtt(std::string date_str, std::string time_str) {
    // convert a UTC date and time string to the corresponding Julian Day of TT

    // extract year, month and day from the date string
    int year, month, day;
    char delim;
    std::stringstream date_ss(date_str);
    date_ss >> day >> delim >> month >> delim >> year;

    // extract hours, minutes and seconds
    int hour, minute, sec;
    std::stringstream time_ss(time_str);
    time_ss >> hour >> delim >> minute >> delim >> sec;

    // construct day fraction
    double day_fraction = (hour / 24.0) + (minute / (24.0 * 60.0)) + (sec / (24.0 * 60.0 * 60.0));

    // convert date to JD UTC
    double utc_int, utc_frac;
    eraCal2jd(year, month, day, &utc_int, &utc_frac);

    // add the day fraction
    utc_frac += day_fraction;
    if (utc_frac >= 1.0) {
        utc_frac -= 1.0;
        utc_int += 1.0;
    }

    // convert UTC to TAI
    double tai_int, tai_frac;
    eraUtctai(utc_int, utc_frac, &tai_int, &tai_frac);

    // convert TAI to TT
    double tt_int, tt_frac;
    eraTaitt(tai_int, tai_frac, &tt_int, &tt_frac);

    // return JD TT
    return tt_int + tt_frac;
}

double utc_date_and_time_to_j2000(std::string date_str, std::string time_str) {
    // convert a UTC dat eand time string to seconds since the J2000.0 epoch

    double jdtt = utc_date_and_time_to_jdtt(date_str, time_str);

    double j2000_seconds = (jdtt - 2451545.0) *86400.0 ; // convert to seconds since J2000.0

    return j2000_seconds;
}



double compute_sun_zenith(double lon, double lat, const Body &earth, const Body &sun) {
    // at lon, lat on earths surface what is the zenith of the sun? // account for oblateness of earth in position vector
    const double& a_earth = earth.Geta();
    const double& b_earth = earth.Getb();

    double a2 = a_earth*a_earth, b2 = b_earth*b_earth;
    double cos_lat = cos(lat), sin_lat = sin(lat);
    double Z = sqrt(a2*cos_lat*cos_lat + b2*sin_lat*sin_lat);

    Vector pos_on_earth = Vector(a2*cos_lat*cos(lon), a2*cos_lat*sin(lon), b2*sin_lat) /Z;
    const Quaternion& q_earth = earth.Getq();
    pos_on_earth = q_earth.rotate(pos_on_earth); // rotate from terrestrial frame to ICRF

    Vector earth_to_sun = sun.Getx() - earth.Getx();
    double cos_zen = pos_on_earth.dot(earth_to_sun/earth_to_sun.norm());

    return M_PI/2 - acos(cos_zen);
}

std::array<double, 3> compute_spherical_seen_from_earth(const Body &body, const Body &earth) {
    // what is the bodies distance, lon and lat as seen from the terrestrial frame?

    Vector r = body.Getx() - earth.Getx();
    const Quaternion& q_earth = earth.Getq();
    r = q_earth.conjugate().rotate(r); // rotate from ICRF to terrestrial frame

    double dist = r.norm();
    double lon = atan2(r(1), r(0));
    double lat = atan(r(2)/sqrt(r(0)*r(0) + r(1)*r(1)));

    return {dist, lon, lat};
}

double compute_shadow_axis_distance_to_earth_center(const Body &earth, const Body &moon, const Body &sun) {
    // the shadow's central axis is given be the line connecting sun and moon
    // return the shortest distance of this line to the earths center

    const Vector& x_earth_ICRF = earth.Getx();
    const Vector& x_moon_ICRF = moon.Getx();
    const Vector& x_sun_ICRF = sun.Getx();
    const Vector& v_moon = moon.Getv();
    const Vector& v_sun = sun.Getv();

    // shift ICRF origin to earth
    // x_earth = 0;
    Vector x_moon = x_moon_ICRF - x_earth_ICRF;
    Vector x_sun = x_sun_ICRF - x_earth_ICRF;

    // we need the retarded positions of the sun and moon (ie corrected for light travel time)
    double PHYS_c = 299792.458; // remove before merge!!!!
    double t_moon = -1.0/PHYS_c *x_moon.norm();
    double t_sun = -1.0/PHYS_c *x_sun.norm();
    x_moon += v_moon*t_moon;
    x_sun += v_sun*t_sun;

    Vector moon_to_sun = x_sun - x_moon;
    moon_to_sun /= moon_to_sun.norm();

    double dist2 = x_moon.squaredNorm() - (moon_to_sun.dot(x_moon))*(moon_to_sun.dot(x_moon));

    return sqrt(dist2);
}

double acos_continued(double x) {
    // extend the acos function by setting it to 0 for x>1 and pi for x>-1
    if (x > 1.0) {return 0.0;}
    else if (x < -1.0) {return M_PI;}
    else {return acos(x);}
}

double atan3(double y, double x) {
    // compute atan2() and cast it to between 0 to 2*pi
    double at2 = atan2(y, x);
    if (at2 >= 0) {return at2;}
    else {return at2 + 2*M_PI;} // if at2 < 0 add 2*pi to make it positive
}

double disks_intersection_area(double R1, double R2, double d) {
    // what is the intersection area of two disks (one with radius R1 and the other with R2)
    // whose centers are separated by a distance d?

    double Z = (R1*R1-R2*R2+d*d)/2/d;
    double angle1 = acos_continued(Z/R1);
    double angle2 = acos_continued((d-Z)/R2);

    return R1*R1*(angle1 - sin(2*angle1)/2) + R2*R2*(angle2 - sin(2*angle2)/2);
}



double cap_solid_angle(double theta) {
    // return the solid angle covered by a spherical cap of angular radius theta
    return 2*M_PI *(1.0-cos(theta));
}


double caps_intersection_solid_angle(double theta1, double theta2, double angle_sep) {
    // what is the intersection solid angle of two spherical caps (one with angular radius theta1 and the other with theta2)
    // whose centers are separated by an angular separation angle_sep?

    if (angle_sep > theta1 + theta2) { // separation is too large: no intersection at all
        return 0.0;
    }
    else if (angle_sep < abs(theta1 - theta2)) { // separation is small enough for total overlap
        return std::min(cap_solid_angle(theta1), cap_solid_angle(theta2));
    }

    // semi-angles of the sectors (phi1 and phi2)
    double cos_phi1 = (cos(theta2) - cos(angle_sep) * cos(theta1)) / (sin(angle_sep) * sin(theta1));
    double cos_phi2 = (cos(theta1) - cos(angle_sep) * cos(theta2)) / (sin(angle_sep) * sin(theta2));
    double phi1 = acos_continued(cos_phi1);
    double phi2 = acos_continued(cos_phi2);

    // spherical Excess using L'Huilier's formula
    double s = (theta1 + theta2 + angle_sep)/2; // Semi-perimeter
    double tan_s2     = tan(s / 2);
    double tan_s_th1  = tan((s - theta1)/2);
    double tan_s_th2  = tan((s - theta2)/2);
    double tan_s_d    = tan((s - angle_sep)/2);
    double tan_E_div_4 = sqrt(tan_s2 * tan_s_th1 * tan_s_th2 * tan_s_d);
    double E = 4*atan(tan_E_div_4);

    // put everything together
    double area1 = 2*phi1 * (1.0 - cos(theta1));
    double area2 = 2*phi2 * (1.0 - cos(theta2));
    double intersection_solid_angle = area1 + area2 - 2*E;

    if (intersection_solid_angle < 0) {
        return 0.0; // if due to numerical precision errors the result is slightly negative return 0
    }

    return intersection_solid_angle;
}



template<size_t grid_size>
void compute_local_occultations(const Body &earth, const Body &moon,  const Body &sun,
                                const heap_array<double, grid_size> &lon_grid, const heap_array<double, grid_size> &lat_grid,
                                heap_array<double, grid_size> &occult_buffer, heap_array<uint8_t, grid_size> &classif_buffer) {
    // compute the occultation ratios of the sun at different points on the earth
    // thereby ignore day and night, i.e. pretend that one can look through the earth and see the eclipse even though its night
    // (this is just so darkness due to eclipse and darkness due to night are not mixed up)
    // write the results into the provided buffers
    const double& a_earth = earth.Geta();
    const double& b_earth = earth.Getb();
    const double& R_moon = moon.GetR();
    const double& R_sun = sun.GetR();

    const Vector& x_earth_ICRF = earth.Getx();
    const Vector& x_moon_ICRF = moon.Getx();
    const Vector& x_sun_ICRF = sun.Getx();
    const Vector& v_moon = moon.Getv();
    const Vector& v_sun = sun.Getv();
    const Quaternion& q_earth = earth.Getq();

    // shift ICRF origin to earth
    // x_earth = 0;
    Vector x_moon = x_moon_ICRF - x_earth_ICRF;
    Vector x_sun = x_sun_ICRF - x_earth_ICRF;

    // rotate into terrestrial frame
    x_moon = q_earth.conjugate().rotate(x_moon);
    x_sun = q_earth.conjugate().rotate(x_sun);

    // we need the retarded positions of the sun and moon (ie corrected for light travel time)
    double PHYS_c = 299792.458; // remove before merge!!!!
    double t_moon = -1.0/PHYS_c *x_moon.norm();
    double t_sun = -1.0/PHYS_c *x_sun.norm();
    x_moon += v_moon*t_moon;
    x_sun += v_sun*t_sun;

    // loop over the longitudes and latitudes and compute the respective occultation rates
    // and the topology including angle of the moon relative to the sun
    #pragma omp parallel for
    for (unsigned i = 0; i < grid_size; i++) {
        double lon = lon_grid[i];
        double lat = lat_grid[i];

        // convert to position vector
        double a2 = a_earth*a_earth, b2 = b_earth*b_earth;
        double cos_lat = cos(lat), sin_lat = sin(lat);
        double Z = sqrt(a2*cos_lat*cos_lat + b2*sin_lat*sin_lat);
        Vector x = Vector(a2*cos_lat*cos(lon), a2*cos_lat*sin(lon), b2*sin_lat) /Z;

        Vector x_s = x_sun - x; // vector from observation point to the sun
        Vector x_m = x_moon - x; // vector from observation point to the moon
        double d_s = x_s.norm(); // distance from observation point to the sun
        double d_m = x_m.norm(); // distance from observation point to the moon
        Vector e_s = x_s/d_s; // unit vector pointing to the sun
        Vector e_m = x_m/d_m; // unit vector pointing to the moon
        double rho_s = sqrt(x_s(0)*x_s(0) + x_s(1)*x_s(1)); // cylindrical radius of x_s
        Vector e_southnorth = Vector(0.0, 0.0, d_s/rho_s) - x_s(2)/rho_s *e_s; // unit vector pointing north to south
        Vector e_eastwest(x_s(1)/rho_s, -x_s(0)/rho_s, 0.0); // unit vector pointing east to west

        double theta_moon = asin(R_moon/d_m); // angular radius of the moon as seen from observation point
        double theta_sun = asin(R_sun/d_s); // angular radius of the sun as seen from observation point

        double angle_sep = acos(e_s.dot(e_m)); // angular separation between the moon and the sun

        // compute intersection solid angle of sun and moon and divide by the sun's solid angle to get occultations
        occult_buffer[i] = caps_intersection_solid_angle(theta_moon, theta_sun, angle_sep) /cap_solid_angle(theta_sun);

        // classify the eclipse topology: no eclipse(0), partial(1), annular(2), total(3)
        uint8_t topology;
        if (angle_sep > theta_sun + theta_moon) {
            topology = 0; // no eclipse
        }
        else if (angle_sep > abs(theta_sun - theta_moon)) {
            topology = 1; // partiaĺ eclipse
        }
        else if (theta_sun > theta_moon) {
            topology = 2; // annular eclipse
        }
        else { // theta_moon >= theta_sun
            topology = 3; // total eclipse
        }

        // calculate the clockwise angle of the moon relative to the sun (0° being north)
        // and sort it into 6° bins
        double d_moon_southnorth = e_southnorth.dot(x_m - x_s);
        double d_moon_eastwest = e_eastwest.dot(x_m - x_s);

        double moon_angle = atan3(d_moon_eastwest, d_moon_southnorth);
        uint8_t moon_angle_binned = floor(moon_angle*30.0/M_PI);

        // combine topology and angle into a combined classification
        classif_buffer[i] = 4*moon_angle_binned + topology;
    }
}



template<size_t grid_size, size_t general_size>
class Eclipse_NetCDF {
public:
    Eclipse_NetCDF() {
        // just instantiate the class, but no file yet
        file_active_ = false;

        // do some sanity checks
        if (general_size == 0) {
            throw std::invalid_argument("NetCDF: general_size must not be 0, since at least time information should be given.");
        }
        if (grid_size == 0) {
            throw std::invalid_argument("NetCDF: grid_size must no be 0. What's the point of all this, if no occultation information is provided?");
        }
    }
    ~Eclipse_NetCDF() {
        // make sure any potentially active file is closed properly before this class is destroyed
        if (ncid_ >= 0) {nc_close(ncid_);}
    }

    bool is_open() const {
        return file_active_;

    }

    void create_new_file(std::string file_path, std::array<std::string, general_size> general_data_keys,
                    const heap_array<double, grid_size> &lon_grid, const heap_array<double, grid_size> &lat_grid) {

        // check for active file
        if (file_active_) {
            throw std::runtime_error("please close the currently active file before creating a new one.");
        }
        // create a new NetCDF file
        bool status = nc_create(file_path.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid_);
        if (status != NC_NOERR) {
            throw std::runtime_error("NetCDF: Failed to create file: " + file_path);
        }

        // define needed dimensions
        nc_def_dim(ncid_, "steps", NC_UNLIMITED, &step_dim_);
        nc_def_dim(ncid_, "grid", grid_size, &grid_dim_);

        // define needed variables
        for (size_t i = 0; i < general_size; i++) {
            nc_def_var(ncid_, general_data_keys[i].c_str(), NC_DOUBLE, 1, &step_dim_, &general_ids_[i]);
        }
        nc_def_var(ncid_, "lon_grid", NC_DOUBLE, 1, &grid_dim_, &lon_grid_id_);
        nc_def_var(ncid_, "lat_grid", NC_DOUBLE, 1, &grid_dim_, &lat_grid_id_);
        int matrix_dim[2] = {step_dim_, grid_dim_};
        nc_def_var(ncid_, "classification_data", NC_UBYTE, 2, matrix_dim, &classif_id_);
        nc_def_var(ncid_, "occultation_data", NC_USHORT, 2, matrix_dim, &occult_id_);
        // add scale_factor and add_offset attributes to occultation_data
        double occult_scale_factor = 1.0 / (pow(2,16)-1);
        double occult_add_offset = 0.0;
        nc_put_att_double(ncid_, occult_id_, "scale_factor", NC_DOUBLE, 1, &occult_scale_factor);
        nc_put_att_double(ncid_, occult_id_, "add_offset", NC_DOUBLE, 1, &occult_add_offset);

        // enable chunking on steps axis
        size_t chunk_sizes[2] = {1, grid_size};
        nc_def_var_chunking(ncid_, classif_id_, NC_CHUNKED, chunk_sizes);
        nc_def_var_chunking(ncid_, occult_id_, NC_CHUNKED, chunk_sizes);

        // enable compression: shuffle=1, deflate=1, level=5
        nc_def_var_deflate(ncid_, classif_id_, 1, 1, 5);
        nc_def_var_deflate(ncid_, occult_id_, 1, 1, 5);

        // close file definition mode
        nc_enddef(ncid_);

        // initialize the start_ variables
        general_start_ = 0;
        matrix_start_[0] = 0, matrix_start_[1] = 0;

        // write the lon and lat grids
        nc_put_var_double(ncid_, lon_grid_id_, lon_grid.data());
        nc_put_var_double(ncid_, lat_grid_id_, lat_grid.data());

        // set file_active flag
        file_active_ = true;
    }

    void close_file() {
        // close file
        if (ncid_ >= 0) {nc_close(ncid_);}

        // reset the file_active flag
        file_active_ = false;
    }

    void write_step(const std::array<double, general_size> &general_data,
                    const heap_array<double, grid_size> &occult_data,
                    const heap_array<uint8_t, grid_size> &classif_data) {

        // check for active file
        if (!file_active_) {
            throw std::runtime_error("NetCDF: cannot write when no active file is present.");
        }

        // write the general data
        for (size_t i = 0; i < general_size; i++) {
            nc_put_vara_double(ncid_, general_ids_[i], &general_start_, &general_count_, &general_data[i]);
        }
        general_start_ += general_count_; // increment the start tracker

        // write the classification data
        nc_put_vara_uchar(ncid_, classif_id_, matrix_start_, matrix_count_, classif_data.data());

        // before writing the occultation data first compress it into the uint16 representation
        for (size_t i = 0; i < grid_size; i++) {
            occult_write_buffer_[i] = static_cast<uint16_t>(occult_data[i] *(pow(2,16)-1) + 0.5);
        }
        // write the compressed occultation data
        nc_put_vara_ushort(ncid_, occult_id_, matrix_start_, matrix_count_, occult_write_buffer_.data());

        matrix_start_[0] += matrix_count_[0]; // increment the start tracker
    }

private:
    int ncid_;
    int step_dim_, grid_dim_;
    std::array<int, general_size> general_ids_;
    size_t general_start_;
    const size_t general_count_{1};
    int lon_grid_id_, lat_grid_id_;
    int occult_id_, classif_id_;
    size_t matrix_start_[2];
    const size_t matrix_count_[2]{1, grid_size};
    heap_array<uint16_t, grid_size> occult_write_buffer_;
    bool file_active_;
};







#endif //UTILS_H_INCLUDED

