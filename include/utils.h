#include <array>
#include <stdexcept>
#include <cmath>
#include <string>

#include <erfa.h>

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED


double PHYS_G = 6.6743015E+04; // km3/(1e24 kg)/s^2, Newtons constant of gravity


class Vector {

public:

    Vector() :
        data_{}
    {}

    Vector(double x, double y, double z) :
        data_{x, y, z}
    {}

    Vector(const Vector &other) :
        data_(other.data_)
    {}

    Vector(std::array<double, 3> data) :
        data_(data)
    {}

    const std::array<double, 3>& GetData() const {
        return data_;
    }

    double& operator[](unsigned idx) {
        return data_[idx];
    }

    const double& operator[](unsigned idx) const {
        return data_[idx];
    }

    void operator=(const Vector &other) {
        data_ = other.data_;
    }

    Vector operator+(const Vector &other) const {
        Vector result;
        for (unsigned i = 0; i < 3; i++) {
            result[i] = data_[i] + other.data_[i];
        }
        return result;
    }

    void operator+=(const Vector &other) {
        *this = *this + other;
    }

    Vector operator-(const Vector &other) const {
        Vector result;
        for (unsigned i = 0; i < 3; i++) {
            result[i] = data_[i] - other.data_[i];
        }
        return result;
    }

    void operator-=(const Vector &other) {
        *this = *this - other;
    }

    double operator*(const Vector &other) const {
        double result = 0.0;
        for (unsigned i = 0; i < 3; i++) {
            result += data_[i] * other.data_[i];
        }
        return result;
    }

    Vector operator*(double scalar) const {
        Vector result;
        for (unsigned i = 0; i < 3; i++) {
            result[i] = data_[i] * scalar;
        }
        return result;
    }

    void operator*=(double scalar) {
        *this = (*this) *scalar;
    }

    Vector operator/(double scalar) const {
        Vector result;
        for (unsigned i = 0; i < 3; i++) {
            result[i] = data_[i] /scalar;
        }
        return result;
    }

    void operator/=(double scalar) {
        *this = (*this) /scalar;
    }

    double Abs() const {
        double result = 0.0;
        for (unsigned i = 0; i < 3; i++) {
            result += data_[i] * data_[i];
        }
        return sqrt(result);
    }

    Vector Hat() const {
        return (*this)/this->Abs();
    }

    Vector CrossProduct(const Vector &other) const {
        Vector result;
        result[0] = data_[1]*other.data_[2] - data_[2]*other.data_[1];
        result[1] = data_[2]*other.data_[0] - data_[0]*other.data_[2];
        result[2] = data_[0]*other.data_[1] - data_[1]*other.data_[0];
        return result;
    }


private:
    std::array<double, 3> data_;
};

Vector operator*(double scalar, const Vector &vec) {
    return vec * scalar;
}


class Matrix {

public:

    Matrix() :
        data_{}
    {}

    Matrix(std::array<double, 3> row_x, std::array<double, 3> row_y, std::array<double, 3> row_z) :
        data_{row_x, row_y, row_z}
    {}

    Matrix(const Matrix &other) :
        data_(other.data_)
    {}

    Matrix(std::array<std::array<double, 3>, 3> data):
        data_(data)
    {}

    const std::array<std::array<double, 3>, 3>& GetData() const {
        return data_;
    }

    std::array<double, 3>& operator[](unsigned row_idx) {
        return data_[row_idx];
    }

    const std::array<double, 3>& operator[](unsigned row_idx) const {
        return data_[row_idx];
    }

    void operator=(const Matrix &other) {
        data_ = other.data_;
    }

    Matrix operator+(const Matrix &other) const {
        Matrix result;
        for (unsigned i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result[i][j] = data_[i][j] + other.data_[i][j];
            }
        }
        return result;
    }

    void operator+=(const Matrix &other) {
        *this = *this + other;
    }

    Matrix operator-(const Matrix &other) const {
        Matrix result;
        for (unsigned i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result[i][j] = data_[i][j] - other.data_[i][j];
            }
        }
        return result;
    }

    void operator-=(const Matrix &other) {
        *this = *this - other;
    }

    Matrix operator*(const Matrix &other) const {
        Matrix result;
        for (unsigned i = 0; i < 3; i++) {
            for (unsigned j = 0; j < 3; j++) {
                for (unsigned k = 0; k < 3; k++) {
                    result[i][j] += data_[i][k] * other.GetData()[k][j];
                }
            }
        }
        return result;
    }

    void operator*=(const Matrix & other) {
        *this = (*this) *other;
    }

    Vector operator*(const Vector &other) const {
        Vector result;
        for (unsigned i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result[i] += data_[i][j] * other.GetData()[j];
            }
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result;
        for (unsigned i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result[i][j] = data_[i][j] * scalar;
            }
        }
        return result;
    }

    void operator*=(double scalar) {
        *this = (*this) *scalar;
    }

    Matrix operator/(double scalar) const {
        Matrix result;
        for (unsigned i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result[i][j] = data_[i][j] / scalar;
            }
        }
        return result;
    }

    void operator/=(double scalar) {
        *this = (*this) /scalar;
    }


private:
    std::array<std::array<double, 3>, 3> data_;
};

Matrix operator*(double scalar, const Matrix &mat) {
    return mat * scalar;
}

Matrix TensorProduct(const Vector &vector_l, const Vector &vector_r) {
    Matrix result;
    for (unsigned i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = vector_l[i] * vector_r[j];
        }
    }
    return result;
}

Matrix Identity{{1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0},
                {0.0, 0.0, 1.0}};


double P2(double x) {
    // second Legendre Polynomial
    return (3*x*x - 1.0)/2;
}



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
        const double& delta = orient_[2];

        Matrix I_inv = {{1.0/Iz_ + 1.0/Ixy_*cos(delta)*cos(delta)/(1.0+sin(delta))/(1.0+sin(delta)), 1.0/Ixy_/(1.0+sin(delta)), 0.0},
                        {1.0/Ixy_/(1.0+sin(delta)), 1.0/Ixy_/cos(delta)/cos(delta), 0.0},
                        {0.0, 0.0, 1.0/Ixy_}};

        return I_inv *L_;
    }


    Vector GetAxis() const {
        const double& alpha = orient_[1];
        const double& delta = orient_[2];

        return Vector(cos(delta)*cos(alpha), cos(delta)*sin(alpha), sin(delta));
    }

    double GetT_rot() const {
        double w = Getw()[0];

        return 2*M_PI /w;
    }

    Matrix GetRotMat() const {
        // get the rotation matrix that transforms from terrestrial coordinates to celestial ones
        const double& phi = orient_[0];
        const double& alpha = orient_[1];
        const double& delta = orient_[2];

        Matrix RotMat_phi = {{cos(phi), -sin(phi), 0.0},
                             {sin(phi), cos(phi), 0.0},
                             {0.0, 0.0, 1.0}};
        Matrix RotMat_axis = {{1 - cos(alpha)*cos(alpha)*(1.0-sin(delta)), -cos(alpha)*sin(alpha)*(1.0-sin(delta)), cos(delta)*cos(alpha)},
                              {-cos(alpha)*sin(alpha)*(1.0-sin(delta)), 1 - sin(alpha)*sin(alpha)*(1.0-sin(delta)), cos(delta)*sin(alpha)},
                              {-cos(delta)*cos(alpha), -cos(delta)*sin(alpha), sin(delta)}};

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
    Vector orient_; // orientation: [0]: rotation angle, [1]: pole RA, [2]: pole Dec
    Vector L_; // angular momentum: components corresponding to orientation vector

    void enforce_orient_range_() {
        double& alpha = orient_[1];
        double& delta = orient_[2];

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



template<typename T, unsigned N>
class Array {
    template<typename T1, unsigned N1>
    friend class Array;
    using TypeArray = Array<T, N>;
    using ScalarArray = Array<double, N>;

public:

    Array() :
        data_{}
    {}

    Array(const TypeArray &other) :
        data_(other.data_)
    {}

    Array(std::array<T, N> data) :
        data_(data)
    {}

    const std::array<T, N>& GetData() const {
        return data_;
    }

    T& operator[](unsigned n) {
        return data_.at(n);
    }

    const T& operator[](unsigned n) const {
        return data_.at(n);
    }

    TypeArray operator+(const TypeArray &other) const {
        TypeArray result;
        for (unsigned i = 0; i < N; i++) {
            result[i] = data_[i] + other.data_[i];
        }
        return result;
    }

    void operator+=(const Array &other) {
        *this = *this + other;
    }

    TypeArray operator-(const TypeArray &other) const {
        TypeArray result;
        for (unsigned i = 0; i < N; i++) {
            result[i] = data_[i] - other.data_[i];
        }
        return result;
    }

    void operator-=(const Array &other) {
        *this = *this - other;
    }

    TypeArray operator*(double scalar) const {
        TypeArray result;
        for (unsigned i = 0; i < N; i++) {
            result[i] = data_[i] * scalar;
        }
        return result;
    }

    void operator*=(double scalar) {
        *this = (*this) *scalar;
    }

    TypeArray operator*(const ScalarArray &scalar_arr) const {
        TypeArray result;
        for (unsigned i = 0; i < N; i++) {
            result[i] = data_[i] * scalar_arr.data_[i];
        }
        return result;
    }

    void operator*=(const ScalarArray &scalar_arr) {
        *this = (*this) *scalar_arr;
    }

    TypeArray operator/(double scalar) const {
        TypeArray result;
        for (unsigned i = 0; i < N; i++) {
            result[i] = data_[i] /scalar;
        }
        return result;
    }

    void operator/=(double scalar) {
        *this = (*this) /scalar;
    }

    TypeArray operator/(const ScalarArray &scalar_arr) const {
        TypeArray result;
        for (unsigned i = 0; i < N; i++) {
            result[i] = data_[i] / scalar_arr.data_[i];
        }
        return result;
    }

    void operator/=(const ScalarArray &scalar_arr) {
        *this = (*this) /scalar_arr;
    }

private:
    std::array<T, N> data_;
};

template<typename T, unsigned N>
Array<T, N> operator*(double scalar, const Array<T, N> &arr) {
    return arr * scalar;
}


template<unsigned N>
using ScalarArray = Array<double, N>;
template<unsigned N>
using VectorArray = Array<Vector, N>;
template<unsigned N>
using MatrixArray = Array<Matrix, N>;


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
            body_arr_[i].Setx(x[i]);
        }
    }

    void Incrementx(const VectorArray<N> &incr_x) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementx(incr_x[i]);
        }
    }

    void Setp(const VectorArray<N> &p) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setp(p[i]);
        }
    }

    void Incrementp(const VectorArray<N> &incr_p) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementp(incr_p[i]);
        }
    }

    void Setorient(const VectorArray<N> &orient) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Setorient(orient[i]);
        }
    }

    void Incrementorient(const VectorArray<N> &incr_orient) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].Incrementorient(incr_orient[i]);
        }
    }

    void SetL(const VectorArray<N> &L) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].SetL(L[i]);
        }
    }

    void IncrementL(const VectorArray<N> &incr_L) {
        for (unsigned i = 0; i < N; i++) {
            body_arr_[i].IncrementL(incr_L[i]);
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

    Vector e_parallel = (x2-x0)/(x2-x0).Abs();

    double d01 = (x1-x0)*e_parallel;

    if (d01 < 0 || d01 > (x2-x0).Abs()) {
        return false;
    }

    double D = (x1-x0 - d01*e_parallel).Abs();
    double r_lightcone_1 = R0 + (R2-R0)*d01/(x2-x0).Abs();

    double A_intersect = intersection_area(R1, r_lightcone_1, D);

    return A_intersect / (M_PI*R0*R0);
}



bool eclipsed(const Vector &x0, double R0, const Vector &x1, double R1, const Vector &x2, double R2) {
    // check if body 1 ecplipses between bodies 0 and 2

    Vector e_parallel = (x2-x0)/(x2-x0).Abs();

    double d01 = (x1-x0)*e_parallel;

    if (d01 < 0 || d01 > (x2-x0).Abs()) {
        return false;
    }

    double D = (x1-x0 - d01*e_parallel).Abs();
    double r_lightcone_1 = R0 + (R2-R0)*d01/(x2-x0).Abs();

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
    double cos_zen = pos_on_earth * earth_to_sun/earth_to_sun.Abs();

    return M_PI/2 - acos(cos_zen);
}






#endif //UTILS_H_INCLUDED

