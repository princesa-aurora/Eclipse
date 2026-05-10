// initial conditions for solar system simulation
// reference frame is ICRF (International Celestial Reference Frame), origin is the Solar System Barycenter (SSB), reference epoch J2000.0
// note that standard ICRF is use, i.e. with the z-axis aligned with earths rotation axis at J2000.0, not the ecliptic plane normal, as is the standard in JPL's Horizons ephemeris
// this is because rotation axes of planets are exclusively to be found in this frame

// time is in Julian Days of TDB (Barycentric Dynamical Time) with initial time at J2000.0 TDB (Terrestrial Time) (January 1st, 2000)
// for the purpore of this simulation the difference between TDB and TT shall be ignored,
// which looses a couple of milliseconds in accuracy, but that's ok
// So we have t0 = t_J2000.0 = 2451545.0 JDTDB =~ 2451545.0 JDTT


#include <include/utils.h>

#ifndef INITIAL_CONDITIONS_H_INCLUDED
#define INITIAL_CONDITIONS_H_INCLUDED


// translational data:
// taken from NASA JPL Horizons System (https://ssd.jpl.nasa.gov/horizons/app.html#/)
// units are: x[km], t[s], v[km/s], M[E+24kg], GM[km^3/s^2]
// R is volumetric mean radius

Vector x0_sun(-1.067706805381631E+06, -3.960361847950951E+05, -1.380651842864990E+05);
Vector v0_sun(9.312571926508239E-03, -1.170150612818866E-02, -5.251266205204772E-03);
double M_sun(1988410.0);
double GM_sun(132712440041.93938);
double R_sun(695700.0);
Vector p0_sun = M_sun * v0_sun;

Vector x0_mercury(-2.052943316392625E+07, -6.032400395765506E+07, -3.013083786350743E+07);
Vector v0_mercury(3.700430442865286E+01, -8.541376791212787E+00, -8.398372410524413E+00);
double M_mercury(0.3302);
double GM_mercury(22031.86855);
double R_mercury(2439.4);
Vector p0_mercury = M_mercury * v0_mercury;

Vector x0_venus(-1.085242008576727E+08, -7.318564957348876E+06, 3.548121862388247E+06);
Vector v0_venus(1.391218600360602E+00, -3.202951993786435E+01, -1.449708673944476E+01);
double M_venus(4.8685);
double GM_venus(324858.592);
double R_venus(6051.84);
Vector p0_venus = M_venus * v0_venus;

Vector x0_earth(-2.756674048064499E+07, 1.323613811539150E+08, 5.741865328641246E+07);
Vector v0_earth(-2.978494749858966E+01, -5.029753814524049E+00, -2.180645068860803E+00);
double M_earth(5.97219);
double GM_earth(398600.435436);
double R_earth(6371.01);
Vector p0_earth = M_earth*v0_earth;

Vector x0_moon(-2.785834886487951E+07, 1.320946643201093E+08, 5.734255079912778E+07);
Vector v0_moon(-2.914141610973326E+01, -5.695841497599683E+00, -2.481970773555618E+00);
double M_moon(0.07349);
double GM_moon(4902.800066);
double R_moon(1737.53);
Vector p0_moon = M_moon * v0_moon;

Vector x0_mars(2.069804338363758E+08, -1.864170129960323E+05, -5.667227498237504E+06);
Vector v0_mars(1.171984975915371E+00, 2.390670819298864E+01, 1.093392065055535E+01);
double M_mars(0.64171);
double GM_mars(42828.375662);
double R_mars(3389.92);
Vector p0_mars = M_mars * v0_mars;

Vector x0_jupiter(5.974999178522581E+08, 4.089902697986262E+08, 1.607562616929507E+08);
Vector v0_jupiter(-7.900547720232828E+00, 1.017187257623615E+01, 4.552504127786899E+00);
double M_jupiter(1898.19);
double GM_jupiter(126686531.900);
double R_jupiter(69911.0);
Vector p0_jupiter = M_jupiter * v0_jupiter;

Vector x0_saturn(9.573176521108806E+08, 9.233194350570333E+08, 3.401627932738827E+08);
Vector v0_saturn(-7.421900386834246E+00, 6.098450820885574E+00, 2.837547973277585E+00);
double M_saturn(568.34);
double GM_saturn(37931206.234);
double R_saturn(58232.0);
Vector p0_saturn = M_saturn * v0_saturn;

Vector x0_uranus(2.157907112723079E+09, -1.871307099571492E+09, -8.501069259962776E+08);
Vector v0_uranus(4.646584677612326E+00, 4.251110198226733E+00, 1.796121552064784E+00);
double M_uranus(86.813);
double GM_uranus(5793950.6103);
double R_uranus(25362.0);
Vector p0_uranus = M_uranus * v0_uranus;

Vector x0_neptune(2.513978764682013E+09, -3.438170142382629E+09, -1.469851583337212E+09);
Vector v0_neptune(4.474587749877780E+00, 2.876585045457400E+00, 1.065773234175670E+00);
double M_neptune(102.409);
double GM_neptune(6835099.97);
double R_neptune(24624.0);
Vector p0_neptune = M_neptune * v0_neptune;


// rotational data::
// taken from IAU Report (https://www.usgs.gov/publications/report-iau-working-group-cartographic-coordinates-and-rotational-elements-2015)
// units are alpha[rad]: pole RA, delta[rad]: pole DEC, W[rad]: prime meridian angle, w[rad/s]: angular velocity
// a[km]: equatorial radius, b[km]: polar radius (taken from IAU 2015, if same as volumetric then Horizons), I[1E+24kg*km^2]: moment of inertia, f[1]: flattness (squares)
// moments of inertia factors taken from (http://astro.vaporia.com/start/momentofinertiafactor.html)

double alpha0_sun(286.13*M_PI/180.0);
double delta0_sun(63.87*M_PI/180.0);
double W0_sun(84.176*M_PI/180.0);
double w0_sun(14.1844000*M_PI/180.0/86400.0);
double a_sun(R_sun); // no noticable oblateness
double b_sun(R_sun);
double I_sun = 0.070 *M_sun*a_sun*a_sun;
double f_sun = 1 - (b_sun*b_sun)/(a_sun*a_sun);

double alpha0_mercury(281.0103*M_PI/180.0);
double delta0_mercury(61.4155*M_PI/180.0);
double W0_mercury(329.5988*M_PI/180.0);
double w0_mercury(6.138510798336938*M_PI/180.0/86400.0);
double a_mercury(2440.53);
double b_mercury(2438.26);
double I_mercury = 0.346 *M_mercury*a_mercury*a_mercury;
double f_mercury = 1 - (b_mercury*b_mercury)/(a_mercury*a_mercury);

double alpha0_venus(272.76*M_PI/180.0);
double delta0_venus(67.16*M_PI/180.0);
double W0_venus(160.20*M_PI/180.0);
double w0_venus(-1.4813688*M_PI/180.0/86400.0);
double a_venus(R_venus); // no noticable oblateness
double b_venus(R_venus);
double I_venus = 0.337 *M_venus*a_venus*a_venus;
double f_venus = 1 - (b_venus*b_venus)/(a_venus*a_venus);

double alpha0_earth(0.0); // by definition of ICRF
double delta0_earth(90.0*M_PI/180.0); // irrelevant since alpha0_earth=0.0
double W0_earth(190.147*M_PI/180.0); // taken from the 2009 IAU report (https://www.usgs.gov/publications/report-iau-working-group-cartographic-coordinates-and-rotational-elements-2009)
double w0_earth(0.00007292115); // taken from Horizons
double a_earth(6378.1366);
double b_earth(6356.7519);
double I_earth = 0.3307 *M_earth*a_earth*a_earth;
double f_earth = 1 - (b_earth*b_earth)/(a_earth*a_earth);

double alpha0_moon(266.9326060656449*M_PI/180.0); // all moon data taken from the 2009 report
double delta0_moon(65.64567138250604*M_PI/180.0);
double W0_moon(41.236695349984274*M_PI/180.0);
double w0_moon(13.178249354527678*M_PI/180.0/86400.0);
double a_moon(R_moon); // no noticable oblateness
double b_moon(R_moon);
double I_moon = 0.3929 *M_moon*a_moon*a_moon;
double f_moon = 1 - (b_moon*b_moon)/(a_moon*a_moon);

double alpha0_mars(317.269202*M_PI/180.0);
double delta0_mars(54.432516*M_PI/180.0);
double W0_mars(176.049863*M_PI/180.0);
double w0_mars(350.89198243554165*M_PI/180.0/86400.0);
double a_mars(3396.19);
double b_mars(3376.20);
double I_mars = 0.3644 *M_mars*a_mars*a_mars;
double f_mars = 1 - (b_mars*b_mars)/(a_mars*a_mars);

double alpha0_jupiter(268.05677758224186*M_PI/180.0);
double delta0_jupiter(64.49489190378664*M_PI/180.0);
double W0_jupiter(284.95*M_PI/180.0);
double w0_jupiter(870.5360000*M_PI/180.0/86400.0);
double a_jupiter(71492.0);
double b_jupiter(66854.0);
double I_jupiter = 0.2756 *M_jupiter*a_jupiter*a_jupiter;
double f_jupiter = 1 - (b_jupiter*b_jupiter)/(a_jupiter*a_jupiter);

double alpha0_saturn(40.589*M_PI/180.0);
double delta0_saturn(83.537*M_PI/180.0);
double W0_saturn(38.90*M_PI/180.0);
double w0_saturn(810.7939024*M_PI/180.0/86400.0);
double a_saturn(60268.0);
double b_saturn(54364.0);
double I_saturn = 0.2234 *M_saturn*a_saturn*a_saturn;
double f_saturn = 1 - (b_saturn*b_saturn)/(a_saturn*a_saturn);

double alpha0_uranus(257.311*M_PI/180.0);
double delta0_uranus(-15.175*M_PI/180.0);
double W0_uranus(203.81*M_PI/180.0);
double w0_uranus(-501.1600928*M_PI/180.0/86400.0);
double a_uranus(25559.0);
double b_uranus(24973.0);
double I_uranus = 0.23 *M_uranus*a_uranus*a_uranus;
double f_uranus = 1 - (b_uranus*b_uranus)/(a_uranus*a_uranus);

double alpha0_neptune(299.33373895877264*M_PI/180.0);
double delta0_neptune(42.95035902184485*M_PI/180.0);
double W0_neptune(249.99600757112736*M_PI/180.0);
double w0_neptune(541.1397756915903*M_PI/180.0/86400.0);
double a_neptune(24764.0);
double b_neptune(24341.0);
double I_neptune = 0.23 *M_neptune*a_neptune*a_neptune;
double f_neptune = 1 - (b_neptune*b_neptune)/(a_neptune*a_neptune);



#endif // INITIAL_CONDITIONS_H_INCLUDED

