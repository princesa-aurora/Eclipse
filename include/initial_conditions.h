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
Vector x0_sun(-1.067706805380953E+06, -3.960361847959462E+05, -1.380651842868809E+05);
Vector v0_sun(9.312571926520472E-03, -1.170150612817771E-02, -5.251266205200356E-03);
double GM_sun(132712440041.279419);
double M_sun = GM_sun /PHYS_G;

double alpha0_sun(286.13*M_PI/180.0);
double delta0_sun(63.87*M_PI/180.0);
double W0_sun(84.176*M_PI/180.0);
double alpha_dot0_sun(0.0);
double delta_dot0_sun(0.0);
double W_dot0_sun(14.1844000*M_PI/180.0/86400.0);
double a1_sun(695700.0);
double a2_sun(695700.0);
double b_sun(695690.0);
double i_f_sun(0.070);
double J2_sun(2.2e-7);
double J22_sun(0.0);
double l2_sun(0.0); // irrelevant due to J22 = 0
std::pair<Quaternion, Quaternion> quaternions_sun = angles_to_quaternions(W0_sun, alpha0_sun, delta0_sun,
                                                                W_dot0_sun, alpha_dot0_sun, delta_dot0_sun);
Quaternion q0_sun = quaternions_sun.first;
Quaternion w0_sun = quaternions_sun.second;

Body Sun("Sun",
        M_sun, a1_sun, a2_sun, b_sun,
        i_f_sun, J2_sun, J22_sun, l2_sun,
        x0_sun, v0_sun,
        q0_sun, w0_sun);

// Mercury
Vector x0_mercury(-2.052943316123468E+07, -6.032400395827633E+07, -3.013083786411830E+07);
Vector v0_mercury(3.700430442920571E+01, -8.541376789510446E+00, -8.398372409672424E+00);
double GM_mercury(22031.868551);
double M_mercury = GM_mercury /PHYS_G;

double alpha0_mercury(281.0103*M_PI/180.0);
double delta0_mercury(61.4155*M_PI/180.0);
double W0_mercury(329.5988*M_PI/180.0);
double alpha_dot0_mercury(-0.0328*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_mercury(-0.0049*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_mercury(6.138510798336938*M_PI/180.0/86400.0);
double a1_mercury(2440.53);
double a2_mercury(2439.28);
double b_mercury(2438.26);
double i_f_mercury(0.346);
double J2_mercury(6.0e-5);
double J22_mercury(1.61e-5);
double l2_mercury(-15.0*M_PI/180.0);
std::pair<Quaternion, Quaternion> quaternions_mercury = angles_to_quaternions(W0_mercury, alpha0_mercury, delta0_mercury,
                                                                W_dot0_mercury, alpha_dot0_mercury, delta_dot0_mercury);
Quaternion q0_mercury = quaternions_mercury.first;
Quaternion w0_mercury = quaternions_mercury.second;

Body Mercury("Mercury",
        M_mercury, a1_mercury, a2_mercury, b_mercury,
        i_f_mercury, J2_mercury, J22_mercury, l2_mercury,
        x0_mercury, v0_mercury,
        q0_mercury, w0_mercury);

// Venus
Vector x0_venus(-1.085242008575715E+08, -7.318564959678600E+06, 3.548121861333776E+06);
Vector v0_venus(1.391218601189967E+00, -3.202951993781091E+01, -1.449708673947320E+01);
double GM_venus(324858.592000);
double M_venus = GM_venus /PHYS_G;

double alpha0_venus(272.76*M_PI/180.0);
double delta0_venus(67.16*M_PI/180.0);
double W0_venus(160.20*M_PI/180.0);
double alpha_dot0_venus(0.0);
double delta_dot0_venus(0.0);
double W_dot0_venus(-1.4813688*M_PI/180.0/86400.0);
double a1_venus(6051.89);
double a2_venus(6051.87);
double b_venus(6051.88);
double i_f_venus(0.337);
double J2_venus(4.4e-6);
double J22_venus(1.01e-6);
double l2_venus(-38.0*M_PI/180.0);
std::pair<Quaternion, Quaternion> quaternions_venus = angles_to_quaternions(W0_venus, alpha0_venus, delta0_venus,
                                                                W_dot0_venus, alpha_dot0_venus, delta_dot0_venus);
Quaternion q0_venus = quaternions_venus.first;
Quaternion w0_venus = quaternions_venus.second;

Body Venus("Venus",
        M_venus, a1_venus, a2_venus, b_venus,
        i_f_venus, J2_venus, J22_venus, l2_venus,
        x0_venus, v0_venus,
        q0_venus, w0_venus);

// Earth (Terra)
Vector x0_earth(-2.756674048281145E+07, 1.323613811535491E+08, 5.741865328625385E+07);
Vector v0_earth(-2.978494749851088E+01, -5.029753814928081E+00, -2.180645069035755E+00);
double GM_earth(398600.435507);
double M_earth = GM_earth /PHYS_G;

double alpha0_earth(0.0); // irrelevant since delta0_earth=90°
double delta0_earth(90.0*M_PI/180.0); // by definition of ICRF
double W0_earth(190.147*M_PI/180.0); // taken from the 2009 IAU report (https://www.usgs.gov/publications/report-iau-working-group-cartographic-coordinates-and-rotational-elements-2009)
double alpha_dot0_earth(-0.641*M_PI/180.0/100/365.25/24/60/60); // taken from the 2009 report
double delta_dot0_earth(-0.557*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_earth(360.986*M_PI/180.0/86400.0); // taken from Horizons
double a1_earth(6378.137);
double a2_earth(6378.1368);
double b_earth(6356.752);
double i_f_earth(0.3307);
double J2_earth(1.08263e-3);
double J22_earth(1.816e-6);
double l2_earth(-14.9*M_PI/180.0);
std::pair<Quaternion, Quaternion> quaternions_earth = angles_to_quaternions(W0_earth, alpha0_earth, delta0_earth,
                                                                W_dot0_earth, alpha_dot0_earth, delta_dot0_earth);
Quaternion q0_earth = quaternions_earth.first;
Quaternion w0_earth = quaternions_earth.second;

Body Earth("Earth",
        M_earth, a1_earth, a2_earth, b_earth,
        i_f_earth, J2_earth, J22_earth, l2_earth,
        x0_earth, v0_earth,
        q0_earth, w0_earth);

// Moon (Luna)
Vector x0_moon(-2.785834886699916E+07, 1.320946643196950E+08, 5.734255079894725E+07);
Vector v0_moon(-2.914141610952193E+01, -5.695841497884601E+00, -2.481970773696848E+00);
double GM_moon(4902.800118);
double M_moon = GM_moon /PHYS_G;

double alpha0_moon(266.9326060656449*M_PI/180.0); // moon orientation data taken from the 2009 report
double delta0_moon(65.64567138250604*M_PI/180.0);
double W0_moon(41.236695349984274*M_PI/180.0);
double alpha_dot0_moon(-78.0062009653193*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_moon(41.1390*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_moon(13.178249354527678*M_PI/180.0/86400.0);
double a1_moon(1738.14);
double a2_moon(1737.33);
double b_moon(1735.97);
double i_f_moon(0.3929);
double J2_moon(2.027e-4);
double J22_moon(1.34e-5);
double l2_moon(0.0);
std::pair<Quaternion, Quaternion> quaternions_moon = angles_to_quaternions(W0_moon, alpha0_moon, delta0_moon,
                                                                W_dot0_moon, alpha_dot0_moon, delta_dot0_moon);
Quaternion q0_moon = quaternions_moon.first;
Quaternion w0_moon = quaternions_moon.second;

Body Moon("Moon",
        M_moon, a1_moon, a2_moon, b_moon,
        i_f_moon, J2_moon, J22_moon, l2_moon,
        x0_moon, v0_moon,
        q0_moon, w0_moon);

// Mars
Vector x0_mars(2.069804338364514E+08, -1.864170114371795E+05, -5.667227497526179E+06);
Vector v0_mars(1.171985008531777E+00, 2.390670819417074E+01, 1.093392063330765E+01);
double GM_mars(42828.375816);
double M_mars = GM_mars /PHYS_G;

double alpha0_mars(317.269202*M_PI/180.0);
double delta0_mars(54.432516*M_PI/180.0);
double W0_mars(176.049863*M_PI/180.0);
double alpha_dot0_mars(-0.10927544735750416*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_mars(-0.058271048430924566*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_mars(350.89198243554165*M_PI/180.0/86400.0);
double a1_mars(3398.0);
double a2_mars(3394.4);
double b_mars(3376.2);
double i_f_mars(0.3644);
double J2_mars(1.960e-3);
double J22_mars(6.31e-5);
double l2_mars(-15.0*M_PI/180.0);
std::pair<Quaternion, Quaternion> quaternions_mars = angles_to_quaternions(W0_mars, alpha0_mars, delta0_mars,
                                                                W_dot0_mars, alpha_dot0_mars, delta_dot0_mars);
Quaternion q0_mars = quaternions_mars.first;
Quaternion w0_mars = quaternions_mars.second;

Body Mars("Mars",
        M_mars, a1_mars, a2_mars, b_mars,
        i_f_mars, J2_mars, J22_mars, l2_mars,
        x0_mars, v0_mars,
        q0_mars, w0_mars);

// Jupiter
Vector x0_jupiter(5.974998767925479E+08, 4.089903139317586E+08, 1.607562819387201E+08);
Vector v0_jupiter(-7.900525116640771E+00, 1.017179630923791E+01, 4.552467787262923E+00);
double GM_jupiter(126712764.100000);
double M_jupiter = GM_jupiter /PHYS_G;

double alpha0_jupiter(268.05677758224186*M_PI/180.0);
double delta0_jupiter(64.49489190378664*M_PI/180.0);
double W0_jupiter(284.95*M_PI/180.0);
double alpha_dot0_jupiter(-0.027569058147491887*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_jupiter(-0.0023646862246096134*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_jupiter(870.5360000*M_PI/180.0/86400.0);
double a1_jupiter(71492.0);
double a2_jupiter(71492.0);
double b_jupiter(66854.0);
double i_f_jupiter(0.2756);
double J2_jupiter(1.4736e-2);
double J22_jupiter(0.0);
double l2_jupiter(0.0); // irrelevent since J22=0
std::pair<Quaternion, Quaternion> quaternions_jupiter = angles_to_quaternions(W0_jupiter, alpha0_jupiter, delta0_jupiter,
                                                                W_dot0_jupiter, alpha_dot0_jupiter, delta_dot0_jupiter);
Quaternion q0_jupiter = quaternions_jupiter.first;
Quaternion w0_jupiter = quaternions_jupiter.second;

Body Jupiter("Jupiter",
        M_jupiter, a1_jupiter, a2_jupiter, b_jupiter,
        i_f_jupiter, J2_jupiter, J22_jupiter, l2_jupiter,
        x0_jupiter, v0_jupiter,
        q0_jupiter, w0_jupiter);

// Saturn
Vector x0_saturn(9.573174174143425E+08, 9.233196218969914E+08, 3.401628003886153E+08);
Vector v0_saturn(-7.422709426014511E+00, 6.097474815228996E+00, 2.837682288255575E+00);
double GM_saturn(37940584.841800);
double M_saturn = GM_saturn /PHYS_G;

double alpha0_saturn(40.589*M_PI/180.0);
double delta0_saturn(83.537*M_PI/180.0);
double W0_saturn(38.90*M_PI/180.0);
double alpha_dot0_saturn(-0.036*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_saturn(-0.004*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_saturn(810.7939024*M_PI/180.0/86400.0);
double a1_saturn(60268.0);
double a2_saturn(60268.0);
double b_saturn(54364.0);
double i_f_saturn(0.2234);
double J2_saturn(1.6291e-2);
double J22_saturn(0.0);
double l2_saturn(0.0); // irrelevant since J22=0
std::pair<Quaternion, Quaternion> quaternions_saturn = angles_to_quaternions(W0_saturn, alpha0_saturn, delta0_saturn,
                                                                W_dot0_saturn, alpha_dot0_saturn, delta_dot0_saturn);
Quaternion q0_saturn = quaternions_saturn.first;
Quaternion w0_saturn = quaternions_saturn.second;

Body Saturn("Saturn",
        M_saturn, a1_saturn, a2_saturn, b_saturn,
        i_f_saturn, J2_saturn, J22_saturn, l2_saturn,
        x0_saturn, v0_saturn,
        q0_saturn, w0_saturn);

// Uranus
Vector x0_uranus(2.157907312953845E+09, -1.871306838939559E+09, -8.501068000312823E+08);
Vector v0_uranus(4.646336807878125E+00, 4.251152675974153E+00, 1.796172785811120E+00);
double GM_uranus(5794556.400000);
double M_uranus = GM_uranus /PHYS_G;

double alpha0_uranus(257.311*M_PI/180.0);
double delta0_uranus(-15.175*M_PI/180.0);
double W0_uranus(203.81*M_PI/180.0);
double alpha_dot0_uranus(0.0);
double delta_dot0_uranus(0.0);
double W_dot0_uranus(-501.1600928*M_PI/180.0/86400.0);
double a1_uranus(25559.0);
double a2_uranus(25559.0);
double b_uranus(24973.0);
double i_f_uranus(0.23);
double J2_uranus(3.341e-3);
double J22_uranus(0.0);
double l2_uranus(0.0); // irrelevant since J22=0
std::pair<Quaternion, Quaternion> quaternions_uranus = angles_to_quaternions(W0_uranus, alpha0_uranus, delta0_uranus,
                                                                W_dot0_uranus, alpha_dot0_uranus, delta_dot0_uranus);
Quaternion q0_uranus = quaternions_uranus.first;
Quaternion w0_uranus = quaternions_uranus.second;

Body Uranus("Uranus",
        M_uranus, a1_uranus, a2_uranus, b_uranus,
        i_f_uranus, J2_uranus, J22_uranus, l2_uranus,
        x0_uranus, v0_uranus,
        q0_uranus, w0_uranus);

// Neptune
Vector x0_neptune(2.513978721723721E+09, -3.438170140316856E+09, -1.469851523010959E+09);
Vector v0_neptune(4.475214621751308E+00, 2.877104855637858E+00, 1.066200548145841E+00);
double GM_neptune(6836527.100580);
double M_neptune = GM_neptune /PHYS_G;

double alpha0_neptune(299.33373895877264*M_PI/180.0);
double delta0_neptune(42.95035902184485*M_PI/180.0);
double W0_neptune(249.99600757112736*M_PI/180.0);
double alpha_dot0_neptune(0.6387097073456971*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_neptune(-0.017470114244133583*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_neptune(541.1397756915903*M_PI/180.0/86400.0);
double a1_neptune(24764.0);
double a2_neptune(24764.0);
double b_neptune(24341.0);
double i_f_neptune(0.23);
double J2_neptune(3.408e-3);
double J22_neptune(0.0);
double l2_neptune(0.0); // irrelevant since J22=0
std::pair<Quaternion, Quaternion> quaternions_neptune = angles_to_quaternions(W0_neptune, alpha0_neptune, delta0_neptune,
                                                                W_dot0_neptune, alpha_dot0_neptune, delta_dot0_neptune);
Quaternion q0_neptune = quaternions_neptune.first;
Quaternion w0_neptune = quaternions_neptune.second;

Body Neptune("Neptune",
        M_neptune, a1_neptune, a2_neptune, b_neptune,
        i_f_neptune, J2_neptune, J22_neptune, l2_neptune,
        x0_neptune, v0_neptune,
        q0_neptune, w0_neptune);

// Pluto
Vector x0_pluto(-1.478399372814445E+09, -4.185975131407770E+09, -8.608780505448110E+08);
Vector v0_pluto(5.253463931238531E+00, -1.964080598078450E+00, -2.195768154117838E+00);
double GM_pluto(975.500000);
double M_pluto = GM_pluto /PHYS_G;

double alpha0_pluto(132.993*M_PI/180.0);
double delta0_pluto(-6.163*M_PI/180.0);
double W0_pluto(302.695*M_PI/180.0);
double alpha_dot0_pluto(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_pluto(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_pluto(-56.3625221*M_PI/180.0/86400.0);
double a1_pluto(1188.7);
double a2_pluto(1188.1);
double b_pluto(1187.1);
double i_f_pluto(0.30);
double J2_pluto(0.0);
double J22_pluto(3.5e-5);
double l2_pluto(0.0);
std::pair<Quaternion, Quaternion> quaternions_pluto = angles_to_quaternions(W0_pluto, alpha0_pluto, delta0_pluto,
                                                                W_dot0_pluto, alpha_dot0_pluto, delta_dot0_pluto);
Quaternion q0_pluto = quaternions_pluto.first;
Quaternion w0_pluto = quaternions_pluto.second;

Body Pluto("Pluto",
        M_pluto, a1_pluto, a2_pluto, b_pluto,
        i_f_pluto, J2_pluto, J22_pluto, l2_pluto,
        x0_pluto, v0_pluto,
        q0_pluto, w0_pluto);

// Eris
Vector x0_eris(1.322247535751656E+10, 5.775062507665586E+09, -1.750958305968350E+09);
Vector v0_eris(-3.431890913481263E-01, 9.396333336817613E-01, 2.047076336369011E+00);
double GM_eris(1108.0);
double M_eris = GM_eris /PHYS_G;

double alpha0_eris(203.4*M_PI/180.0);
double delta0_eris(-61.6*M_PI/180.0);
double W0_eris(0.0*M_PI/180.0);
double alpha_dot0_eris(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_eris(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_eris(22.80424*M_PI/180.0/86400.0);
double a1_eris(1163.0);
double a2_eris(1163.0);
double b_eris(1163.0);
double i_f_eris(0.370);
double J2_eris(0.0);
double J22_eris(0.0);
double l2_eris(0.0); // irrelevant since J22=0
std::pair<Quaternion, Quaternion> quaternions_eris = angles_to_quaternions(W0_eris, alpha0_eris, delta0_eris,
                                                                W_dot0_eris, alpha_dot0_eris, delta_dot0_eris);
Quaternion q0_eris = quaternions_eris.first;
Quaternion w0_eris = quaternions_eris.second;

Body Eris("Eris",
        M_eris, a1_eris, a2_eris, b_eris,
        i_f_eris, J2_eris, J22_eris, l2_eris,
        x0_eris, v0_eris,
        q0_eris, w0_eris);

// Ceres
Vector x0_ceres(-3.570100537503446E+08, 8.123520469843985E+07, 1.107476694834988E+08);
Vector v0_ceres(-6.196624728846403E+00, -1.705827086161680E+01, -6.765800523726184E+00);
double GM_ceres(62.62890);
double M_ceres = GM_ceres /PHYS_G;

double alpha0_ceres(291.4274*M_PI/180.0);
double delta0_ceres(66.7603*M_PI/180.0);
double W0_ceres(170.309*M_PI/180.0);
double alpha_dot0_ceres(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_ceres(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_ceres(952.1532635*M_PI/180.0/86400.0);
double a1_ceres(483.1);
double a2_ceres(481.0);
double b_ceres(445.9);
double i_f_ceres(0.37);
double J2_ceres(0.02583);
double J22_ceres(1.81e-3);
double l2_ceres(0.0);
std::pair<Quaternion, Quaternion> quaternions_ceres = angles_to_quaternions(W0_ceres, alpha0_ceres, delta0_ceres,
                                                                W_dot0_ceres, alpha_dot0_ceres, delta_dot0_ceres);
Quaternion q0_ceres = quaternions_ceres.first;
Quaternion w0_ceres = quaternions_ceres.second;

Body Ceres("Ceres",
        M_ceres, a1_ceres, a2_ceres, b_ceres,
        i_f_ceres, J2_ceres, J22_ceres, l2_ceres,
        x0_ceres, v0_ceres,
        q0_ceres, w0_ceres);

// Vesta
Vector x0_vesta(-2.035604403464110E+08, -2.428277492720578E+08, -7.020476895727536E+07);
Vector v0_vesta(1.667424655301909E+01, -1.112020471863626E+01, -6.606096271993724E+00);
double GM_vesta(17.288245);
double M_vesta = GM_vesta /PHYS_G;

double alpha0_vesta(309.031*M_PI/180.0);
double delta0_vesta(42.235*M_PI/180.0);
double W0_vesta(285.39*M_PI/180.0);
double alpha_dot0_vesta(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_vesta(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_vesta(1617.3329428*M_PI/180.0/86400.0);
double a1_vesta(286.3);
double a2_vesta(278.6);
double b_vesta(223.2);
double i_f_vesta(0.320);
double J2_vesta(0.03178);
double J22_vesta(3.0e-3);
double l2_vesta(0.0);
std::pair<Quaternion, Quaternion> quaternions_vesta = angles_to_quaternions(W0_vesta, alpha0_vesta, delta0_vesta,
                                                                W_dot0_vesta, alpha_dot0_vesta, delta_dot0_vesta);
Quaternion q0_vesta = quaternions_vesta.first;
Quaternion w0_vesta = quaternions_vesta.second;

Body Vesta("Vesta",
        M_vesta, a1_vesta, a2_vesta, b_vesta,
        i_f_vesta, J2_vesta, J22_vesta, l2_vesta,
        x0_vesta, v0_vesta,
        q0_vesta, w0_vesta);

// Pallas
Vector x0_pallas(-1.269002268927843E+08, 2.904887908103527E+08, -4.912457343501391E+07);
Vector v0_pallas(-2.031372264343589E+01, -9.189004891206704E+00, 3.220148192842943E+00);
double GM_pallas(13.63);
double M_pallas = GM_pallas /PHYS_G;

double alpha0_pallas(38.3213*M_PI/180.0);
double delta0_pallas(66.5392*M_PI/180.0);
double W0_pallas(0.0*M_PI/180.0);
double alpha_dot0_pallas(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_pallas(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_pallas(1105.772591*M_PI/180.0/86400.0);
double a1_pallas(275.0);
double a2_pallas(258.0);
double b_pallas(238.0);
double i_f_pallas(0.375);
double J2_pallas(1.51e-2);
double J22_pallas(1.2e-2);
double l2_pallas(0.0);
std::pair<Quaternion, Quaternion> quaternions_pallas = angles_to_quaternions(W0_pallas, alpha0_pallas, delta0_pallas,
                                                                W_dot0_pallas, alpha_dot0_pallas, delta_dot0_pallas);
Quaternion q0_pallas = quaternions_pallas.first;
Quaternion w0_pallas = quaternions_pallas.second;

Body Pallas("Pallas",
        M_pallas, a1_pallas, a2_pallas, b_pallas,
        i_f_pallas, J2_pallas, J22_pallas, l2_pallas,
        x0_pallas, v0_pallas,
        q0_pallas, w0_pallas);

// Hygiea
Vector x0_hygiea(-3.562223996255096E+08, -1.906738360989999E+08, -1.116844471544101E+08);
Vector v0_hygiea(1.058803870658594E+01, -1.440527063339443E+01, -5.761049138115527E+00);
double GM_hygiea(5.78);
double M_hygiea = GM_hygiea /PHYS_G;

double alpha0_hygiea(319*M_PI/180.0);
double delta0_hygiea(-46*M_PI/180.0);
double W0_hygiea(0.0*M_PI/180.0);
double alpha_dot0_hygiea(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_hygiea(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_hygiea(624.93175*M_PI/180.0/86400.0);
double a1_hygiea(217.0);
double a2_hygiea(213.0);
double b_hygiea(197.0);
double i_f_hygiea(0.380);
double J2_hygiea(0.013);
double J22_hygiea(3.0e-3);
double l2_hygiea(0.0);
std::pair<Quaternion, Quaternion> quaternions_hygiea = angles_to_quaternions(W0_hygiea, alpha0_hygiea, delta0_hygiea,
                                                                W_dot0_hygiea, alpha_dot0_hygiea, delta_dot0_hygiea);
Quaternion q0_hygiea = quaternions_hygiea.first;
Quaternion w0_hygiea = quaternions_hygiea.second;

Body Hygiea("Hygiea",
        M_hygiea, a1_hygiea, a2_hygiea, b_hygiea,
        i_f_hygiea, J2_hygiea, J22_hygiea, l2_hygiea,
        x0_hygiea, v0_hygiea,
        q0_hygiea, w0_hygiea);

// Iris
Vector x0_iris(-1.762246893448783E+08, 2.558763868953050E+08, 8.778781118816169E+07);
Vector v0_iris(-1.973641676997331E+01, -6.292921870466768E+00, -4.636583904954921E+00);
double GM_iris(0.9);
double M_iris = GM_iris /PHYS_G;

double alpha0_iris(19.0*M_PI/180.0);
double delta0_iris(19.0*M_PI/180.0);
double W0_iris(0.0*M_PI/180.0);
double alpha_dot0_iris(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_iris(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_iris(50.43*M_PI/180.0/86400.0);
double a1_iris(126.0);
double a2_iris(110.0);
double b_iris(90.0);
double i_f_iris(0.4);
double J2_iris(0.1);
double J22_iris(1.8e-2);
double l2_iris(0.0);
std::pair<Quaternion, Quaternion> quaternions_iris = angles_to_quaternions(W0_iris, alpha0_iris, delta0_iris,
                                                                W_dot0_iris, alpha_dot0_iris, delta_dot0_iris);
Quaternion q0_iris = quaternions_iris.first;
Quaternion w0_iris = quaternions_iris.second;

Body Iris("Iris",
        M_iris, a1_iris, a2_iris, b_iris,
        i_f_iris, J2_iris, J22_iris, l2_iris,
        x0_iris, v0_iris,
        q0_iris, w0_iris);

// Bamberga
Vector x0_bamberga(7.905177457679231E+07, -2.953300496802147E+08, -1.811485559864970E+08);
Vector v0_bamberga(1.727915740023365E+01, 8.046544160781494E+00, 7.310258317029986E+00);
double GM_bamberga(0.61);
double M_bamberga = GM_bamberga /PHYS_G;

double alpha0_bamberga(326.0*M_PI/180.0);
double delta0_bamberga(-37.0*M_PI/180.0);
double W0_bamberga(0.0*M_PI/180.0);
double alpha_dot0_bamberga(0.0*M_PI/180.0/100/365.25/24/60/60);
double delta_dot0_bamberga(0.0*M_PI/180.0/100/365.25/24/60/60);
double W_dot0_bamberga(12.23*M_PI/180.0/86400.0);
double a1_bamberga(119.0);
double a2_bamberga(113.0);
double b_bamberga(102.0);
double i_f_bamberga(0.4);
double J2_bamberga(0.0);
double J22_bamberga(6.0e-3);
double l2_bamberga(0.0);
std::pair<Quaternion, Quaternion> quaternions_bamberga = angles_to_quaternions(W0_bamberga, alpha0_bamberga, delta0_bamberga,
                                                                W_dot0_bamberga, alpha_dot0_bamberga, delta_dot0_bamberga);
Quaternion q0_bamberga = quaternions_bamberga.first;
Quaternion w0_bamberga = quaternions_bamberga.second;

Body Bamberga("Bamberga",
        M_bamberga, a1_bamberga, a2_bamberga, b_bamberga,
        i_f_bamberga, J2_bamberga, J22_bamberga, l2_bamberga,
        x0_bamberga, v0_bamberga,
        q0_bamberga, w0_bamberga);




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
                                        {"Eris", Eris},
                                        {"Ceres", Ceres},
                                        {"Vesta", Vesta},
                                        {"Pallas", Pallas},
                                        {"Hygiea", Hygiea},
                                        {"Iris", Iris},
                                        {"Bamberga", Bamberga}
};



#endif // INITIAL_CONDITIONS_H_INCLUDED

