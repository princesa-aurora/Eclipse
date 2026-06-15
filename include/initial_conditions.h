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
// units are: x[km], t[s], v[km/s], GM[km^3/s^2], M[E+24kg]

// rotational data::
// taken from IAU Report (https://www.usgs.gov/publications/report-iau-working-group-cartographic-coordinates-and-rotational-elements-2015)
// units are alpha[rad]: pole RA, delta[rad]: pole DEC, W[rad]: prime meridian angle, w[rad/s]: angular velocity
// a[km]: equatorial radius, b[km]: polar radius (taken from IAU 2015, if same as volumetric then Horizons), I[1E+24kg*km^2]: moment of inertia, f[1]: flattness (squares)
// moments of inertia factors taken from (http://astro.vaporia.com/start/momentofinertiafactor.html)

// Sun (Sol)
Vector x0_sun(-1.067706805381631E+06, -3.960361847950951E+05, -1.380651842864990E+05);
Vector v0_sun(9.312571926508239E-03, -1.170150612818866E-02, -5.251266205204772E-03);
double GM_sun(132712440041.93938);
double M_sun = GM_sun /PHYS_G;

double alpha0_sun(286.13*M_PI/180.0);
double delta0_sun(63.87*M_PI/180.0);
double W0_sun(84.176*M_PI/180.0);
double phi0_sun = W0_sun + alpha0_sun - M_PI/2;
double alpha_dot0_sun(0.0);
double delta_dot0_sun(0.0);
double W_dot0_sun(14.1844000*M_PI/180.0/86400.0);
double phi_dot0_sun = W_dot0_sun + alpha_dot0_sun;
double a_sun(695700.0); // no noticable oblateness
double b_sun(695700.0);
double i_f_sun = 0.070;
double Iz_sun = i_f_sun *M_sun*a_sun*a_sun;
double Ixy_sun = i_f_sun *M_sun*(a_sun*a_sun + b_sun*b_sun)/2;
double J2_sun = 2.2e-7;
Vector orient0_sun(phi0_sun, alpha0_sun, delta0_sun);
Vector w0_sun(phi_dot0_sun, alpha_dot0_sun, delta_dot0_sun);

Body Sun("Sun",
        M_sun, a_sun, b_sun, i_f_sun, J2_sun,
        x0_sun, v0_sun,
        orient0_sun, w0_sun);

// Mercury
Vector x0_mercury(-2.052943316392625E+07, -6.032400395765506E+07, -3.013083786350743E+07);
Vector v0_mercury(3.700430442865286E+01, -8.541376791212787E+00, -8.398372410524413E+00);
double GM_mercury(22031.86855);
double M_mercury = GM_mercury /PHYS_G;

double alpha0_mercury(281.0103*M_PI/180.0);
double delta0_mercury(61.4155*M_PI/180.0);
double W0_mercury(329.5988*M_PI/180.0);
double phi0_mercury = W0_mercury + alpha0_mercury - M_PI/2;
double alpha_dot0_mercury(-0.0328*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_mercury(-0.0049*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_mercury(6.138510798336938*M_PI/180.0/86400.0);
double phi_dot0_mercury = W_dot0_mercury + alpha_dot0_mercury;
double a_mercury(2439.7);
double b_mercury(2439.7); // no noticable oblateness
double i_f_mercury = 0.346;
double Iz_mercury = i_f_mercury *M_mercury*a_mercury*a_mercury;
double Ixy_mercury = i_f_mercury *M_mercury*(a_mercury*a_mercury + b_mercury*b_mercury)/2;
double J2_mercury = 6.0e-5;
Vector orient0_mercury(phi0_mercury, alpha0_mercury, delta0_mercury);
Vector w0_mercury(phi_dot0_mercury, alpha_dot0_mercury, delta_dot0_mercury);

Body Mercury("Mercury",
        M_mercury, a_mercury, b_mercury, i_f_mercury, J2_mercury,
        x0_mercury, v0_mercury,
        orient0_mercury, w0_mercury);

// Venus
Vector x0_venus(-1.085242008576727E+08, -7.318564957348876E+06, 3.548121862388247E+06);
Vector v0_venus(1.391218600360602E+00, -3.202951993786435E+01, -1.449708673944476E+01);
double GM_venus(324858.592);
double M_venus = GM_venus /PHYS_G;

double alpha0_venus(272.76*M_PI/180.0);
double delta0_venus(67.16*M_PI/180.0);
double W0_venus(160.20*M_PI/180.0);
double phi0_venus = W0_venus + alpha0_venus - M_PI/2;
double alpha_dot0_venus(0.0);
double delta_dot0_venus(0.0);
double W_dot0_venus(-1.4813688*M_PI/180.0/86400.0);
double phi_dot0_venus = W_dot0_venus + alpha_dot0_venus;
double a_venus(6051.84); // no noticable oblateness
double b_venus(6051.84);
double i_f_venus = 0.337;
double Iz_venus = i_f_venus *M_venus*a_venus*a_venus;
double Ixy_venus = i_f_venus *M_venus*(a_venus*a_venus + b_venus*b_venus)/2;
double J2_venus = 4.4e-6;
Vector orient0_venus(phi0_venus, alpha0_venus, delta0_venus);
Vector w0_venus(phi_dot0_venus, alpha_dot0_venus, delta_dot0_venus);

Body Venus("Venus",
        M_venus, a_venus, b_venus, i_f_venus, J2_venus,
        x0_venus, v0_venus,
        orient0_venus, w0_venus);

// Earth (Terra)
Vector x0_earth(-2.756674048064499E+07, 1.323613811539150E+08, 5.741865328641246E+07);
Vector v0_earth(-2.978494749858966E+01, -5.029753814524049E+00, -2.180645068860803E+00);
double GM_earth(398600.435436);
double M_earth = GM_earth /PHYS_G;

double alpha0_earth(0.0); // irrelevant since delta0_earth=90°
double delta0_earth(90.0*M_PI/180.0); // by definition of ICRF
double W0_earth(190.147*M_PI/180.0); // taken from the 2009 IAU report (https://www.usgs.gov/publications/report-iau-working-group-cartographic-coordinates-and-rotational-elements-2009)
double phi0_earth = W0_earth - alpha0_earth + M_PI/2;
double alpha_dot0_earth(-0.641*M_PI/180.0/100/365.25/24/60/60); // taken from the 2009 report
double delta_dot0_earth(-0.557*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_earth(0.00007292115); // taken from Horizons
double phi_dot0_earth = W_dot0_earth + alpha_dot0_earth;
double a_earth(6378.1366);
double b_earth(6356.7519);
double i_f_earth = 0.3307;
double Iz_earth = i_f_earth *M_earth*a_earth*a_earth;
double Ixy_earth = i_f_earth *M_earth*(a_earth*a_earth + b_earth*b_earth)/2;
double J2_earth = 1.08263e-3;
Vector orient0_earth(phi0_earth, alpha0_earth, delta0_earth);
Vector w0_earth(phi_dot0_earth, alpha_dot0_earth, delta_dot0_earth);

Body Earth("Earth",
        M_earth, a_earth, b_earth, i_f_earth, J2_earth,
        x0_earth, v0_earth,
        orient0_earth, w0_earth);

// Moon (Luna)
Vector x0_moon(-2.785834886487951E+07, 1.320946643201093E+08, 5.734255079912778E+07);
Vector v0_moon(-2.914141610973326E+01, -5.695841497599683E+00, -2.481970773555618E+00);
double GM_moon(4902.800066);
double M_moon = GM_moon /PHYS_G;

double alpha0_moon(266.9326060656449*M_PI/180.0); // all moon data taken from the 2009 report
double delta0_moon(65.64567138250604*M_PI/180.0);
double W0_moon(41.236695349984274*M_PI/180.0);
double phi0_moon = W0_moon + alpha0_moon - M_PI/2;
double alpha_dot0_moon(-78.0062009653193*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_moon(41.1390*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_moon(13.178249354527678*M_PI/180.0/86400.0);
double phi_dot0_moon = W_dot0_moon + alpha_dot0_moon;
double a_moon(1738.14);
double b_moon(1735.97);
double i_f_moon = 0.3929;
double Iz_moon = i_f_moon *M_moon*a_moon*a_moon;
double Ixy_moon = i_f_moon *M_moon*(a_moon*a_moon + b_moon*b_moon)/2;
double J2_moon = 2.027e-4;
Vector orient0_moon(phi0_moon, alpha0_moon, delta0_moon);
Vector w0_moon(phi_dot0_moon, alpha_dot0_moon, delta_dot0_moon);

Body Moon("Moon",
        M_moon, a_moon, b_moon, i_f_moon, J2_moon,
        x0_moon, v0_moon,
        orient0_moon, w0_moon);

// Mars
Vector x0_mars(2.069804338363758E+08, -1.864170129960323E+05, -5.667227498237504E+06);
Vector v0_mars(1.171984975915371E+00, 2.390670819298864E+01, 1.093392065055535E+01);
double GM_mars(42828.375662);
double M_mars = GM_mars /PHYS_G;

double alpha0_mars(317.269202*M_PI/180.0);
double delta0_mars(54.432516*M_PI/180.0);
double W0_mars(176.049863*M_PI/180.0);
double phi0_mars = W0_mars + alpha0_mars - M_PI/2;
double alpha_dot0_mars(-0.10927544735750416*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_mars(-0.058271048430924566*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_mars(350.89198243554165*M_PI/180.0/86400.0);
double phi_dot0_mars = W_dot0_mars + alpha_dot0_mars;
double a_mars(3396.19);
double b_mars(3376.20);
double i_f_mars = 0.3644;
double Iz_mars = i_f_mars *M_mars*a_mars*a_mars;
double Ixy_mars = i_f_mars *M_mars*(a_mars*a_mars + b_mars*b_mars)/2;
double J2_mars = 1.960e-3;
Vector orient0_mars(phi0_mars, alpha0_mars, delta0_mars);
Vector w0_mars(phi_dot0_mars, alpha_dot0_mars, delta_dot0_mars);

Body Mars("Mars",
        M_mars, a_mars, b_mars, i_f_mars, J2_mars,
        x0_mars, v0_mars,
        orient0_mars, w0_mars);

// Jupiter
Vector x0_jupiter(5.974998767931225e+08, 4.089903139310188e+08, 1.607562819383890e+08);
Vector v0_jupiter(-7.900525116626677, 10.17179630924754, 4.552467787266711);
double GM_jupiter(126712764.1);
double M_jupiter = GM_jupiter /PHYS_G;

double alpha0_jupiter(268.05677758224186*M_PI/180.0);
double delta0_jupiter(64.49489190378664*M_PI/180.0);
double W0_jupiter(284.95*M_PI/180.0);
double phi0_jupiter = W0_jupiter + alpha0_jupiter - M_PI/2;
double alpha_dot0_jupiter(-0.027569058147491887*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_jupiter(-0.0023646862246096134*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_jupiter(870.5360000*M_PI/180.0/86400.0);
double phi_dot0_jupiter = W_dot0_jupiter + alpha_dot0_jupiter;
double a_jupiter(71492.0);
double b_jupiter(66854.0);
double i_f_jupiter = 0.2756;
double Iz_jupiter = i_f_jupiter *M_jupiter*a_jupiter*a_jupiter;
double Ixy_jupiter = i_f_jupiter *M_jupiter*(a_jupiter*a_jupiter + b_jupiter*b_jupiter)/2;
double J2_jupiter = 1.4736e-2;
Vector orient0_jupiter(phi0_jupiter, alpha0_jupiter, delta0_jupiter);
Vector w0_jupiter(phi_dot0_jupiter, alpha_dot0_jupiter, delta_dot0_jupiter);

Body Jupiter("Jupiter",
        M_jupiter, a_jupiter, b_jupiter, i_f_jupiter, J2_jupiter,
        x0_jupiter, v0_jupiter,
        orient0_jupiter, w0_jupiter);

// Saturn
Vector x0_saturn(9.573174174148824e+08, 9.233196218965478e+08, 3.401628003884089e+08);
Vector v0_saturn(-7.422709426010933, 6.097474815232451, 2.837682288256848);
double GM_saturn(37940584.8418);
double M_saturn = GM_saturn /PHYS_G;

double alpha0_saturn(40.589*M_PI/180.0);
double delta0_saturn(83.537*M_PI/180.0);
double W0_saturn(38.90*M_PI/180.0);
double phi0_saturn = W0_saturn + alpha0_saturn - M_PI/2;
double alpha_dot0_saturn(-0.036*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_saturn(-0.004*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_saturn(810.7939024*M_PI/180.0/86400.0);
double phi_dot0_saturn = W_dot0_saturn + alpha_dot0_saturn;
double a_saturn(60268.0);
double b_saturn(54364.0);
double i_f_saturn = 0.2234;
double Iz_saturn = i_f_saturn *M_saturn*a_saturn*a_saturn;
double Ixy_saturn = i_f_saturn *M_saturn*(a_saturn*a_saturn + b_saturn*b_saturn)/2;
double J2_saturn = 1.6291e-2;
Vector orient0_saturn(phi0_saturn, alpha0_saturn, delta0_saturn);
Vector w0_saturn(phi_dot0_saturn, alpha_dot0_saturn, delta_dot0_saturn);

Body Saturn("Saturn",
        M_saturn, a_saturn, b_saturn, i_f_saturn, J2_saturn,
        x0_saturn, v0_saturn,
        orient0_saturn, w0_saturn);

// Uranus
Vector x0_uranus(2.157907312953506e+09, -1.871306838939868e+09, -8.501068000314130e+08);
Vector v0_uranus(4.646336807878913, 4.251152675973469, 1.796172785810810);
double GM_uranus(5794565.6);
double M_uranus = GM_uranus /PHYS_G;

double alpha0_uranus(257.311*M_PI/180.0);
double delta0_uranus(-15.175*M_PI/180.0);
double W0_uranus(203.81*M_PI/180.0);
double phi0_uranus = W0_uranus + alpha0_uranus - M_PI/2;
double alpha_dot0_uranus(0.0);
double delta_dot0_uranus(0.0);
double W_dot0_uranus(-501.1600928*M_PI/180.0/86400.0);
double phi_dot0_uranus = W_dot0_uranus + alpha_dot0_uranus;
double a_uranus(25559.0);
double b_uranus(24973.0);
double i_f_uranus = 0.23;
double Iz_uranus = i_f_uranus *M_uranus*a_uranus*a_uranus;
double Ixy_uranus = i_f_uranus *M_uranus*(a_uranus*a_uranus + b_uranus*b_uranus)/2;
double J2_uranus = 3.341e-3;
Vector orient0_uranus(phi0_uranus, alpha0_uranus, delta0_uranus);
Vector w0_uranus(phi_dot0_uranus, alpha_dot0_uranus, delta_dot0_uranus);

Body Uranus("Uranus",
        M_uranus, a_uranus, b_uranus, i_f_uranus, J2_uranus,
        x0_uranus, v0_uranus,
        orient0_uranus, w0_uranus);

// Neptune
Vector x0_neptune(2.513978721723395e+09, -3.438170140317066e+09, -1.469851523011037e+09);
Vector v0_neptune(4.475214621751574, 2.877104855637495, 1.066200548145686);
double GM_neptune(6836525.2);
double M_neptune = GM_neptune /PHYS_G;

double alpha0_neptune(299.33373895877264*M_PI/180.0);
double delta0_neptune(42.95035902184485*M_PI/180.0);
double W0_neptune(249.99600757112736*M_PI/180.0);
double phi0_neptune = W0_neptune + alpha0_neptune - M_PI/2;
double alpha_dot0_neptune(0.6387097073456971*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_neptune(-0.017470114244133583*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_neptune(541.1397756915903*M_PI/180.0/86400.0);
double phi_dot0_neptune = W_dot0_neptune + alpha_dot0_neptune;
double a_neptune(24764.0);
double b_neptune(24341.0);
double i_f_neptune = 0.23;
double Iz_neptune = i_f_neptune *M_neptune*a_neptune*a_neptune;
double Ixy_neptune = i_f_neptune *M_neptune*(a_neptune*a_neptune + b_neptune*b_neptune)/2;
double J2_neptune = 3.408e-3;
Vector orient0_neptune(phi0_neptune, alpha0_neptune, delta0_neptune);
Vector w0_neptune(phi_dot0_neptune, alpha_dot0_neptune, delta_dot0_neptune);

Body Neptune("Neptune",
        M_neptune, a_neptune, b_neptune, i_f_neptune, J2_neptune,
        x0_neptune, v0_neptune,
        orient0_neptune, w0_neptune);



#endif // INITIAL_CONDITIONS_H_INCLUDED

