// initial conditions for solar system simulation
// units x[km], t[s], v[km/s], M[E+24kg], GM[km^3/s^2]
// coordinate system center: Solar System Barycenter (SSB) [500@0]
// data taken from NASA JPL Horizons System


#include <include/utils.h>

#ifndef INITIAL_CONDITIONS_H_INCLUDED
#define INITIAL_CONDITIONS_H_INCLUDED


double t0_mjdtdb(0.0); // Modified Julian Day, Barycentric Dynamical Time


Vector x0_sun(-2.147366862601276E+05, -1.107314740516006E+06, 7.872687406769139E+03);
Vector v0_sun(1.497230954044547E-02, -2.286760594509669E-03, -3.852850787786648E-04);
double M_sun(1988410.0);
double GM_sun(132712440041.93938);
double R_sun(695700.0);
Vector p0_sun = M_sun * v0_sun;

Vector x0_mercury(-1.641671016174383E+06, -7.030982145832370E+07, -5.500901691279922E+06);
Vector v0_mercury(3.895712519047996E+01, 1.533540782024240E+00, -3.463667519399485E+00);
double M_mercury(0.3302);
double GM_mercury(22031.86855);
double R_mercury(2439.4);
Vector p0_mercury = M_mercury * v0_mercury;

Vector x0_venus(8.235656966950777E+07, 6.874648584757969E+07, -3.840334655031845E+06);
Vector v0_venus(-2.271278398589050E+01, 2.658328505096891E+01, 1.667017426726334E+00);
double M_venus(4.8685);
double GM_venus(324858.592);
double R_venus(6051.84);
Vector p0_venus = M_venus * v0_venus;

Vector x0_earth(8.167478050099756E+07, 1.219761784301356E+08, 4.943201072294265E+04);
Vector v0_earth(-2.527200819085060E+01, 1.637895052963432E+01, 3.180445933331555E-03);
double M_earth(5.97219);
double GM_earth(398600.435436);
double R_earth(6371.01);
Vector p0_earth = M_earth*v0_earth;

Vector x0_moon(8.205720272882979E+07, 1.219936549911627E+08, 6.391529419735819E+04);
Vector v0_moon(-2.539156844562578E+01, 1.739417870909824E+01, 8.189674945994341E-02);
double M_moon(0.07349);
double GM_moon(4902.800066);
double R_moon(1737.53);
Vector p0_moon = M_moon * v0_moon;

Vector x0_mars(2.008772641065271E+08, -4.993145330508292E+07, -6.016639488632187E+06);
Vector v0_mars(6.683521268023283E+00, 2.559590709738389E+01, 3.685880963570831E-01);
double M_mars(0.64171);
double GM_mars(42828.375662);
double R_mars(3389.92);
Vector p0_mars = M_mars * v0_mars;

Vector x0_jupiter(1.759511361776932E+08, 7.377613014380314E+08, -6.921718972510695E+06);
Vector v0_jupiter(-1.285891899124471E+01, 3.641192846378708E+00, 2.739715763027799E-01);
double M_jupiter(1898.19);
double GM_jupiter(126686531.900);
double R_jupiter(69911.0);
Vector p0_jupiter = M_jupiter * v0_jupiter;

Vector x0_saturn(-8.417178318316665E+08, 1.069674433866184E+09, 1.437907121638811E+07);
Vector v0_saturn(-8.123465164370547E+00, -5.995319469562134E+00, 4.282929236172841E-01);
double M_saturn(568.34);
double GM_saturn(37931206.234);
double R_saturn(58232.0);
Vector p0_saturn = M_saturn * v0_saturn;

Vector x0_uranus(1.304340077732821E+09, 2.595641971154305E+09, -7.255972440584779E+06);
Vector v0_uranus(-6.138582961281127E+00, 2.740180600217014E+00, 9.017622119480506E-02);
double M_uranus(86.813);
double GM_uranus(5793950.6103);
double R_uranus(25362.0);
Vector p0_uranus = M_uranus * v0_uranus;

Vector x0_neptune(4.460344330760460E+09, -3.249820332645947E+08, -9.605593621564491E+07);
Vector v0_neptune(3.583721593800415E-01, 5.452894203339737E+00, -1.202267813343165E-01);
double M_neptune(102.409);
double GM_neptune(6835099.97);
double R_neptune(24624.0);
Vector p0_neptune = M_neptune * v0_neptune;



#endif // INITIAL_CONDITIONS_H_INCLUDED

