// initial conditions for solar system simulation
// reference frame is ICRF (International Celestial Reference Frame), origin is the Solar System Barycenter (SSB), reference epoch J2000.0
// note that standard ICRF is use, i.e. with the z-axis aligned with earths rotation axis at J2000.0, not the ecliptic plane normal, as is the standard in JPL's Horizons ephemeris
// this is because rotation axes of planets are exclusively to be found in this frame

// time is in Julian Days of TDB (Barycentric Dynamical Time) with initial time at J2000.0 TDB (Terrestrial Time) (January 1st, 2000)
// for the purpore of this simulation the difference between TDB and TT shall be ignored,
// which looses a couple of milliseconds in accuracy, but that's ok
// So we have t0 = t_J2000.0 = 2451545.0 JDTDB =~ 2451545.0 JDTT


#include <include/utils.h>
#include <include/physics.h>

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
double GM_sun(132712440041.279419);
double M_sun = GM_sun /PHYS_G;

double alpha0_sun(286.13*M_PI/180.0);
double delta0_sun(63.87*M_PI/180.0);
double W0_sun(84.176*M_PI/180.0);
double alpha_dot0_sun(0.0);
double delta_dot0_sun(0.0);
double W_dot0_sun(14.1844000*M_PI/180.0/86400.0);
double a_sun(695700.0); // no noticable oblateness
double b_sun(695700.0);
double i_f_sun = 0.070;
double J2_sun = 2.2e-7;
std::pair<Quaternion, Quaternion> quaternions_sun = angles_to_quaternions(W0_sun, alpha0_sun, delta0_sun,
                                                                W_dot0_sun, alpha_dot0_sun, delta_dot0_sun);
Quaternion q0_sun = quaternions_sun.first;
Quaternion w0_sun = quaternions_sun.second;

Body Sun("Sun",
        M_sun, a_sun, b_sun, i_f_sun, J2_sun,
        x0_sun, v0_sun,
        q0_sun, w0_sun);

// Mercury
Vector x0_mercury(-2.052943316392625E+07, -6.032400395765506E+07, -3.013083786350743E+07);
Vector v0_mercury(3.700430442865286E+01, -8.541376791212787E+00, -8.398372410524413E+00);
double GM_mercury(22031.868551);
double M_mercury = GM_mercury /PHYS_G;

double alpha0_mercury(281.0103*M_PI/180.0);
double delta0_mercury(61.4155*M_PI/180.0);
double W0_mercury(329.5988*M_PI/180.0);
double alpha_dot0_mercury(-0.0328*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_mercury(-0.0049*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_mercury(6.138510798336938*M_PI/180.0/86400.0);
double a_mercury(2439.7);
double b_mercury(2439.7); // no noticable oblateness
double i_f_mercury = 0.346;
double J2_mercury = 6.0e-5;
std::pair<Quaternion, Quaternion> quaternions_mercury = angles_to_quaternions(W0_mercury, alpha0_mercury, delta0_mercury,
                                                                W_dot0_mercury, alpha_dot0_mercury, delta_dot0_mercury);
Quaternion q0_mercury = quaternions_mercury.first;
Quaternion w0_mercury = quaternions_mercury.second;

Body Mercury("Mercury",
        M_mercury, a_mercury, b_mercury, i_f_mercury, J2_mercury,
        x0_mercury, v0_mercury,
        q0_mercury, w0_mercury);

// Venus
Vector x0_venus(-1.085242008576727E+08, -7.318564957348876E+06, 3.548121862388247E+06);
Vector v0_venus(1.391218600360602E+00, -3.202951993786435E+01, -1.449708673944476E+01);
double GM_venus(324858.592000);
double M_venus = GM_venus /PHYS_G;

double alpha0_venus(272.76*M_PI/180.0);
double delta0_venus(67.16*M_PI/180.0);
double W0_venus(160.20*M_PI/180.0);
double alpha_dot0_venus(0.0);
double delta_dot0_venus(0.0);
double W_dot0_venus(-1.4813688*M_PI/180.0/86400.0);
double a_venus(6051.84); // no noticable oblateness
double b_venus(6051.84);
double i_f_venus = 0.337;
double J2_venus = 4.4e-6;
std::pair<Quaternion, Quaternion> quaternions_venus = angles_to_quaternions(W0_venus, alpha0_venus, delta0_venus,
                                                                W_dot0_venus, alpha_dot0_venus, delta_dot0_venus);
Quaternion q0_venus = quaternions_venus.first;
Quaternion w0_venus = quaternions_venus.second;

Body Venus("Venus",
        M_venus, a_venus, b_venus, i_f_venus, J2_venus,
        x0_venus, v0_venus,
        q0_venus, w0_venus);

// Earth (Terra)
Vector x0_earth(-2.756674048064499E+07, 1.323613811539150E+08, 5.741865328641246E+07);
Vector v0_earth(-2.978494749858966E+01, -5.029753814524049E+00, -2.180645068860803E+00);
double GM_earth(398600.435507);
double M_earth = GM_earth /PHYS_G;

double alpha0_earth(0.0); // irrelevant since delta0_earth=90°
double delta0_earth(90.0*M_PI/180.0); // by definition of ICRF
double W0_earth(190.147*M_PI/180.0); // taken from the 2009 IAU report (https://www.usgs.gov/publications/report-iau-working-group-cartographic-coordinates-and-rotational-elements-2009)
double alpha_dot0_earth(-0.641*M_PI/180.0/100/365.25/24/60/60); // taken from the 2009 report
double delta_dot0_earth(-0.557*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_earth(360.986*M_PI/180.0/86400.0); // taken from Horizons
double a_earth(6378.1366);
double b_earth(6356.7519);
double i_f_earth = 0.3307;
double J2_earth = 1.08263e-3;
std::pair<Quaternion, Quaternion> quaternions_earth = angles_to_quaternions(W0_earth, alpha0_earth, delta0_earth,
                                                                W_dot0_earth, alpha_dot0_earth, delta_dot0_earth);
Quaternion q0_earth = quaternions_earth.first;
Quaternion w0_earth = quaternions_earth.second;

Body Earth("Earth",
        M_earth, a_earth, b_earth, i_f_earth, J2_earth,
        x0_earth, v0_earth,
        q0_earth, w0_earth);

// Moon (Luna)
Vector x0_moon(-2.785834886487951E+07, 1.320946643201093E+08, 5.734255079912778E+07);
Vector v0_moon(-2.914141610973326E+01, -5.695841497599683E+00, -2.481970773555618E+00);
double GM_moon(4902.800118);
double M_moon = GM_moon /PHYS_G;

double alpha0_moon(266.9326060656449*M_PI/180.0); // all moon data taken from the 2009 report
double delta0_moon(65.64567138250604*M_PI/180.0);
double W0_moon(41.236695349984274*M_PI/180.0);
double alpha_dot0_moon(-78.0062009653193*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_moon(41.1390*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_moon(13.178249354527678*M_PI/180.0/86400.0);
double a_moon(1738.14);
double b_moon(1735.97);
double i_f_moon = 0.3929;
double J2_moon = 2.027e-4;
std::pair<Quaternion, Quaternion> quaternions_moon = angles_to_quaternions(W0_moon, alpha0_moon, delta0_moon,
                                                                W_dot0_moon, alpha_dot0_moon, delta_dot0_moon);
Quaternion q0_moon = quaternions_moon.first;
Quaternion w0_moon = quaternions_moon.second;

Body Moon("Moon",
        M_moon, a_moon, b_moon, i_f_moon, J2_moon,
        x0_moon, v0_moon,
        q0_moon, w0_moon);

// Mars
Vector x0_mars(2.069804338363758E+08, -1.864170129960323E+05, -5.667227498237504E+06);
Vector v0_mars(1.171984975915371E+00, 2.390670819298864E+01, 1.093392065055535E+01);
double GM_mars(42828.375816);
double M_mars = GM_mars /PHYS_G;

double alpha0_mars(317.269202*M_PI/180.0);
double delta0_mars(54.432516*M_PI/180.0);
double W0_mars(176.049863*M_PI/180.0);
double alpha_dot0_mars(-0.10927544735750416*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_mars(-0.058271048430924566*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_mars(350.89198243554165*M_PI/180.0/86400.0);
double a_mars(3396.19);
double b_mars(3376.20);
double i_f_mars = 0.3644;
double J2_mars = 1.960e-3;
std::pair<Quaternion, Quaternion> quaternions_mars = angles_to_quaternions(W0_mars, alpha0_mars, delta0_mars,
                                                                W_dot0_mars, alpha_dot0_mars, delta_dot0_mars);
Quaternion q0_mars = quaternions_mars.first;
Quaternion w0_mars = quaternions_mars.second;

Body Mars("Mars",
        M_mars, a_mars, b_mars, i_f_mars, J2_mars,
        x0_mars, v0_mars,
        q0_mars, w0_mars);

// Jupiter
Vector x0_jupiter(5.974998767931225e+08, 4.089903139310188e+08, 1.607562819383890e+08);
Vector v0_jupiter(-7.900525116626677, 10.17179630924754, 4.552467787266711);
double GM_jupiter(126712764.100000);
double M_jupiter = GM_jupiter /PHYS_G;

double alpha0_jupiter(268.05677758224186*M_PI/180.0);
double delta0_jupiter(64.49489190378664*M_PI/180.0);
double W0_jupiter(284.95*M_PI/180.0);
double alpha_dot0_jupiter(-0.027569058147491887*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_jupiter(-0.0023646862246096134*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_jupiter(870.5360000*M_PI/180.0/86400.0);
double a_jupiter(71492.0);
double b_jupiter(66854.0);
double i_f_jupiter = 0.2756;
double J2_jupiter = 1.4736e-2;
std::pair<Quaternion, Quaternion> quaternions_jupiter = angles_to_quaternions(W0_jupiter, alpha0_jupiter, delta0_jupiter,
                                                                W_dot0_jupiter, alpha_dot0_jupiter, delta_dot0_jupiter);
Quaternion q0_jupiter = quaternions_jupiter.first;
Quaternion w0_jupiter = quaternions_jupiter.second;

Body Jupiter("Jupiter",
        M_jupiter, a_jupiter, b_jupiter, i_f_jupiter, J2_jupiter,
        x0_jupiter, v0_jupiter,
        q0_jupiter, w0_jupiter);

// Saturn
Vector x0_saturn(9.573174174148824e+08, 9.233196218965478e+08, 3.401628003884089e+08);
Vector v0_saturn(-7.422709426010933, 6.097474815232451, 2.837682288256848);
double GM_saturn(37940584.841800);
double M_saturn = GM_saturn /PHYS_G;

double alpha0_saturn(40.589*M_PI/180.0);
double delta0_saturn(83.537*M_PI/180.0);
double W0_saturn(38.90*M_PI/180.0);
double alpha_dot0_saturn(-0.036*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_saturn(-0.004*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_saturn(810.7939024*M_PI/180.0/86400.0);
double a_saturn(60268.0);
double b_saturn(54364.0);
double i_f_saturn = 0.2234;
double J2_saturn = 1.6291e-2;
std::pair<Quaternion, Quaternion> quaternions_saturn = angles_to_quaternions(W0_saturn, alpha0_saturn, delta0_saturn,
                                                                W_dot0_saturn, alpha_dot0_saturn, delta_dot0_saturn);
Quaternion q0_saturn = quaternions_saturn.first;
Quaternion w0_saturn = quaternions_saturn.second;

Body Saturn("Saturn",
        M_saturn, a_saturn, b_saturn, i_f_saturn, J2_saturn,
        x0_saturn, v0_saturn,
        q0_saturn, w0_saturn);

// Uranus
Vector x0_uranus(2.157907312953506e+09, -1.871306838939868e+09, -8.501068000314130e+08);
Vector v0_uranus(4.646336807878913, 4.251152675973469, 1.796172785810810);
double GM_uranus(5794556.400000);
double M_uranus = GM_uranus /PHYS_G;

double alpha0_uranus(257.311*M_PI/180.0);
double delta0_uranus(-15.175*M_PI/180.0);
double W0_uranus(203.81*M_PI/180.0);
double alpha_dot0_uranus(0.0);
double delta_dot0_uranus(0.0);
double W_dot0_uranus(-501.1600928*M_PI/180.0/86400.0);
double a_uranus(25559.0);
double b_uranus(24973.0);
double i_f_uranus = 0.23;
double J2_uranus = 3.341e-3;
std::pair<Quaternion, Quaternion> quaternions_uranus = angles_to_quaternions(W0_uranus, alpha0_uranus, delta0_uranus,
                                                                W_dot0_uranus, alpha_dot0_uranus, delta_dot0_uranus);
Quaternion q0_uranus = quaternions_uranus.first;
Quaternion w0_uranus = quaternions_uranus.second;

Body Uranus("Uranus",
        M_uranus, a_uranus, b_uranus, i_f_uranus, J2_uranus,
        x0_uranus, v0_uranus,
        q0_uranus, w0_uranus);

// Neptune
Vector x0_neptune(2.513978721723395e+09, -3.438170140317066e+09, -1.469851523011037e+09);
Vector v0_neptune(4.475214621751574, 2.877104855637495, 1.066200548145686);
double GM_neptune(6836527.100580);
double M_neptune = GM_neptune /PHYS_G;

double alpha0_neptune(299.33373895877264*M_PI/180.0);
double delta0_neptune(42.95035902184485*M_PI/180.0);
double W0_neptune(249.99600757112736*M_PI/180.0);
double alpha_dot0_neptune(0.6387097073456971*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_neptune(-0.017470114244133583*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_neptune(541.1397756915903*M_PI/180.0/86400.0);
double a_neptune(24764.0);
double b_neptune(24341.0);
double i_f_neptune = 0.23;
double J2_neptune = 3.408e-3;
std::pair<Quaternion, Quaternion> quaternions_neptune = angles_to_quaternions(W0_neptune, alpha0_neptune, delta0_neptune,
                                                                W_dot0_neptune, alpha_dot0_neptune, delta_dot0_neptune);
Quaternion q0_neptune = quaternions_neptune.first;
Quaternion w0_neptune = quaternions_neptune.second;

Body Neptune("Neptune",
        M_neptune, a_neptune, b_neptune, i_f_neptune, J2_neptune,
        x0_neptune, v0_neptune,
        q0_neptune, w0_neptune);

// Pluto
Vector x0_pluto(-1.478399372814827E+09, -4.185975131407627E+09, -8.608780505446513E+08);
Vector v0_pluto(5.253463931238377E+00, -1.964080598078888E+00, -2.195768154117928E+00);
double GM_pluto(975.500000);
double M_pluto = GM_pluto /PHYS_G;

double alpha0_pluto(132.993*M_PI/180.0);
double delta0_pluto(-6.163*M_PI/180.0);
double W0_pluto(302.695*M_PI/180.0);
double alpha_dot0_pluto(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_pluto(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_pluto(-56.3625221*M_PI/180.0/86400.0);
double a_pluto(1188.3);
double b_pluto(1188.3);
double i_f_pluto = 0.30;
double J2_pluto = 0.0;
std::pair<Quaternion, Quaternion> quaternions_pluto = angles_to_quaternions(W0_pluto, alpha0_pluto, delta0_pluto,
                                                                W_dot0_pluto, alpha_dot0_pluto, delta_dot0_pluto);
Quaternion q0_pluto = quaternions_pluto.first;
Quaternion w0_pluto = quaternions_pluto.second;

Body Pluto("Pluto",
        M_pluto, a_pluto, b_pluto, i_f_pluto, J2_pluto,
        x0_pluto, v0_pluto,
        q0_pluto, w0_pluto);

// Eris
Vector x0_eris(1.322247535751659E+10, 5.775062507665518E+09, -1.750958305968499E+09);
Vector v0_eris(-3.431890913480846E-01, 9.396333336817794E-01, 2.047076336369006E+00);
double GM_eris(1108.0);
double M_eris = GM_eris /PHYS_G;

double alpha0_eris(203.4*M_PI/180.0);
double delta0_eris(-61.6*M_PI/180.0);
double W0_eris(0.0*M_PI/180.0);
double alpha_dot0_eris(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_eris(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_eris(22.80424*M_PI/180.0/86400.0);
double a_eris(1163.0);
double b_eris(1163.0);
double i_f_eris = 0.370;
double J2_eris = 0.0;
std::pair<Quaternion, Quaternion> quaternions_eris = angles_to_quaternions(W0_eris, alpha0_eris, delta0_eris,
                                                                W_dot0_eris, alpha_dot0_eris, delta_dot0_eris);
Quaternion q0_eris = quaternions_eris.first;
Quaternion w0_eris = quaternions_eris.second;

Body Eris("Eris",
        M_eris, a_eris, b_eris, i_f_eris, J2_eris,
        x0_eris, v0_eris,
        q0_eris, w0_eris);

// Ceres
Vector x0_ceres(-3.570100537498939E+08, 8.123520469968063E+07, 1.107476694839909E+08);
Vector v0_ceres(-6.196624728908223E+00, -1.705827086160263E+01, -6.765800523706930E+00);
double GM_ceres(62.62890);
double M_ceres = GM_ceres /PHYS_G;

double alpha0_ceres(291.4274*M_PI/180.0);
double delta0_ceres(66.7603*M_PI/180.0);
double W0_ceres(170.309*M_PI/180.0);
double alpha_dot0_ceres(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_ceres(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_ceres(952.1532635*M_PI/180.0/86400.0);
double a_ceres(482.0);
double b_ceres(445.9);
double i_f_ceres = 0.37;
double J2_ceres = 0.02583;
std::pair<Quaternion, Quaternion> quaternions_ceres = angles_to_quaternions(W0_ceres, alpha0_ceres, delta0_ceres,
                                                                W_dot0_ceres, alpha_dot0_ceres, delta_dot0_ceres);
Quaternion q0_ceres = quaternions_ceres.first;
Quaternion w0_ceres = quaternions_ceres.second;

Body Ceres("Ceres",
        M_ceres, a_ceres, b_ceres, i_f_ceres, J2_ceres,
        x0_ceres, v0_ceres,
        q0_ceres, w0_ceres);

// Vesta
Vector x0_vesta(-2.035604403476238E+08, -2.428277492712490E+08, -7.020476895679487E+07);
Vector v0_vesta(1.667424655296137E+01, -1.112020471870536E+01, -6.606096272013695E+00);
double GM_vesta(17.288245);
double M_vesta = GM_vesta /PHYS_G;

double alpha0_vesta(309.031*M_PI/180.0);
double delta0_vesta(42.235*M_PI/180.0);
double W0_vesta(285.39*M_PI/180.0);
double alpha_dot0_vesta(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_vesta(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_vesta(1617.3329428*M_PI/180.0/86400.0);
double a_vesta(283.0);
double b_vesta(223.2);
double i_f_vesta = 0.320;
double J2_vesta = 0.03178;
std::pair<Quaternion, Quaternion> quaternions_vesta = angles_to_quaternions(W0_vesta, alpha0_vesta, delta0_vesta,
                                                                W_dot0_vesta, alpha_dot0_vesta, delta_dot0_vesta);
Quaternion q0_vesta = quaternions_vesta.first;
Quaternion w0_vesta = quaternions_vesta.second;

Body Vesta("Vesta",
        M_vesta, a_vesta, b_vesta, i_f_vesta, J2_vesta,
        x0_vesta, v0_vesta,
        q0_vesta, w0_vesta);

// Pallas
Vector x0_pallas(-1.269002268913077E+08, 2.904887908110207E+08, -4.912457343524797E+07);
Vector v0_pallas(-2.031372264347271E+01, -9.189004891121632E+00, 3.220148192828610E+00);
double GM_pallas(13.66);
double M_pallas = GM_pallas /PHYS_G;

double alpha0_pallas(38.3213*M_PI/180.0);
double delta0_pallas(66.5392*M_PI/180.0);
double W0_pallas(0.0*M_PI/180.0);
double alpha_dot0_pallas(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_pallas(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_pallas(1105.772591*M_PI/180.0/86400.0);
double a_pallas(269.0);
double b_pallas(249.0);
double i_f_pallas = 0.375;
double J2_pallas = 0.0151;
std::pair<Quaternion, Quaternion> quaternions_pallas = angles_to_quaternions(W0_pallas, alpha0_pallas, delta0_pallas,
                                                                W_dot0_pallas, alpha_dot0_pallas, delta_dot0_pallas);
Quaternion q0_pallas = quaternions_pallas.first;
Quaternion w0_pallas = quaternions_pallas.second;

Body Pallas("Pallas",
        M_pallas, a_pallas, b_pallas, i_f_pallas, J2_pallas,
        x0_pallas, v0_pallas,
        q0_pallas, w0_pallas);

// Hygiea
Vector x0_hygiea(-3.562223996262794E+08, -1.906738360979529E+08, -1.116844471539913E+08);
Vector v0_hygiea(1.058803870653905E+01, -1.440527063341956E+01, -5.761049138130254E+00);
double GM_hygiea(5.83);
double M_hygiea = GM_hygiea /PHYS_G;

double alpha0_hygiea(319*M_PI/180.0);
double delta0_hygiea(-46*M_PI/180.0);
double W0_hygiea(0.0*M_PI/180.0);
double alpha_dot0_hygiea(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_hygiea(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_hygiea(624.93175*M_PI/180.0/86400.0);
double a_hygiea(220.0);
double b_hygiea(212.0);
double i_f_hygiea = 0.380;
double J2_hygiea = 0.013;
std::pair<Quaternion, Quaternion> quaternions_hygiea = angles_to_quaternions(W0_hygiea, alpha0_hygiea, delta0_hygiea,
                                                                W_dot0_hygiea, alpha_dot0_hygiea, delta_dot0_hygiea);
Quaternion q0_hygiea = quaternions_hygiea.first;
Quaternion w0_hygiea = quaternions_hygiea.second;

Body Hygiea("Hygiea",
        M_hygiea, a_hygiea, b_hygiea, i_f_hygiea, J2_hygiea,
        x0_hygiea, v0_hygiea,
        q0_hygiea, w0_hygiea);




std::map<std::string, Body> initial_bodies_map{
                                        {"Sun", Sun},
                                        {"Mercury", Mercury},
                                        {"Venus", Venus},
                                        {"Earth", Earth},
                                        {"Moon", Moon},
                                        {"Mars", Mars},
                                        {"Jupiter", Jupiter},
                                        {"Saturn", Saturn},
                                        {"Uranus", Uranus},
                                        {"Neptune", Neptune},
                                        {"Pluto", Pluto},
                                        {"Ceres", Ceres},
                                        {"Vesta", Vesta},
                                        {"Pallas", Pallas},
                                        {"Hygiea", Hygiea}
};






#endif // INITIAL_CONDITIONS_H_INCLUDED

