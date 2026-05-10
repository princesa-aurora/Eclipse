#include <array>
#include <stdexcept>
#include <cmath>
#include <string>

#include <erfa.h>

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

class Vector {

public:

    Vector() :
        data_()
    {}

    Vector(double x, double y, double z) :
        data_{x, y, z}
    {}

    Vector(const Vector &other) :
        data_(other.GetData())
    {}

    Vector(std::array<double, 3> data) :
        data_(data)
    {}

    std::array<double, 3> GetData() const {
        return data_;
    }

    double GetX() const {
        return data_[0];
    }

    double GetY() const {
        return data_[1];
    }

    double GetZ() const {
        return data_[2];
    }

    void operator=(const Vector &other) {
        data_ = other.data_;
    }

    Vector operator+(const Vector &other) const {
        std::array<double, 3> result;
        for (int i = 0; i < 3; i++) {
            result[i] = data_[i] + other.data_[i];
        }
        return Vector(result);
    }

    Vector operator-(const Vector &other) const {
        std::array<double, 3> result;
        for (int i = 0; i < 3; i++) {
            result[i] = data_[i] - other.data_[i];
        }
        return Vector(result);
    }

    double operator*(const Vector &other) const {
        double result = 0.0;
        for (int i = 0; i < 3; i++) {
            result += data_[i] * other.data_[i];
        }
        return result;
    }

    Vector operator*(double scalar) const {
        std::array<double, 3> result;
        for (int i = 0; i < 3; i++) {
            result[i] = data_[i] * scalar;
        }
        return Vector(result);
    }

    Vector operator/(double scalar) const {
        std::array<double, 3> result;
        for (int i = 0; i < 3; i++) {
            result[i] = data_[i] /scalar;
        }
        return Vector(result);
    }

    double Abs() const {
        double result = 0.0;
        for (int i = 0; i < 3; i++) {
            result += data_[i] * data_[i];
        }
        return sqrt(result);
    }

private:
    std::array<double, 3> data_;
};

Vector operator*(double scalar, const Vector &vec) {
    return vec * scalar;
}




template<typename T, int N>
class Array {
    template<typename T1, int N1>
    friend class Array;
    using TypeArray = Array<T, N>;
    using ScalarArray = Array<double, N>;

public:

    Array() :
        data_()
    {}

    Array(const TypeArray &other) :
        data_(other.GetData())
    {}

    Array(std::array<T, N> data) :
        data_(data)
    {}

    std::array<T, N> GetData() const {
        return data_;
    }

    T& operator[](int n) {
        return data_.at(n);
    }

    T at(int n) const {
        return data_.at(n);
    }

    TypeArray operator+(const TypeArray &other) const {
        std::array<T, N> result;
        for (int i = 0; i < N; i++) {
            result[i] = data_[i] + other.data_[i];
        }
        return TypeArray(result);
    }

    TypeArray operator-(const TypeArray &other) const {
        std::array<T, N> result;
        for (int i = 0; i < N; i++) {
            result[i] = data_[i] - other.data_[i];
        }
        return TypeArray(result);
    }

    TypeArray operator*(double scalar) const {
        std::array<T, N> result;
        for (int i = 0; i < N; i++) {
            result[i] = data_[i] * scalar;
        }
        return TypeArray(result);
    }

    TypeArray operator*(const ScalarArray &scalar_arr) const {
        std::array<T, N> result;
        for (int i = 0; i < N; i++) {
            result[i] = data_[i] * scalar_arr.data_[i];
        }
        return TypeArray(result);
    }

    TypeArray operator/(double scalar) const {
        std::array<T, N> result;
        for (int i = 0; i < N; i++) {
            result[i] = data_[i] /scalar;
        }
        return TypeArray(result);
    }

    TypeArray operator/(const ScalarArray &scalar_arr) const {
        std::array<T, N> result;
        for (int i = 0; i < N; i++) {
            result[i] = data_[i] / scalar_arr.data_[i];
        }
        return TypeArray(result);
    }

private:
    std::array<T, N> data_;
};

template<typename T, int N>
Array<T, N> operator*(double scalar, const Array<T, N> &arr) {
    return arr * scalar;
}



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


std::string mjdtdb_to_utc_date(double mjdtdb) {
    // convert Modified Julian Date (days) in TDB to a UTC date string

    double tdb = mjdtdb + 2400000.5; // convert MJD to JD
    double tdb_int = floor(tdb); // split into integer
    double tdb_frac = tdb - tdb_int; // and fractional parts

    // convert TDB to TT (Terrestrial Time)
    double jdtt;
    double tt_int, tt_frac;
    eraTdbtt(tdb_int, tdb_frac, 0.0, &tt_int, &tt_frac); // dtr is actually not 0.0 but up to 3.3ms, which we shall ignore here

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
    return std::to_string(day) + "/" + std::to_string(month) + "/" + std::to_string(year) + " " +
           std::to_string(hour) + ":" + std::to_string(minute) + ":" + std::to_string(sec);
}





#endif //UTILS_H_INCLUDED

