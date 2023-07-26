#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#define imax 400


const long double Pi = acos(-1.0);

long double xobs, yobs;
long double epsi3, a13, a22, a52;
long double spin;
long double iobs_deg;
int phicount;

/*---Also set computataional parameters of polish doughnut (Shafqat)*/
/* If I want to change the spin of the BH, then I have to change these four parameters below */
long double wcusp; /* The value of effective potential calculated corresponding to the cusp surface. This is calculated by mathematica code named "combined.nb" */
long double angm_disk;     /* Angular momentum per unit energy of the disk which remains constant (because of considering constant angular momentum disk). also calculated by "combined.nb"*/
long double disk_length; /* length of the disk. Also calculated from "combined.nb"*/ 
long double disk_inner_edge; /* the radial distance from BH to the position of cusp also calculated by "combined.nb" */
/*-----------------------------------------------------------*/


void raytrace(long double xobs, long double yobs, long double iobs, long double xin, long double disk_length_combined ,long double traces[], int& stop_integration);
void diffeqs(long double b, long double vars[], long double diffs[]);
void redshift(long double r, long double ktkp, long double& gg);
//void redshift_polish_doughnut(long double r, long double th, long double l ,long double ktkp, long double& gg);
void intersection(long double x_1, long double y_1, long double z_1, long double x_2, long double y_2, long double z_2, long double x_d[]);
void metric(long double z1, long double z2, long double mn[][4]);
void metric_rderivatives(long double z1, long double z2, long double dmn[][4]);
void find_isco(long double z1, long double& isco);
//void polish_doughnut(long double r, long double theta, long double phi ,long double angm_disk, long double spin, long double& w_current);
//void emission_angle_rth(long double r, long double th, long double l, long double a, long double kr, long double kphi ,long double kth , long double lambda , long double& em_angle);

#include "diffeqs.cpp"
#include "intersection.cpp"
#include "metric.cpp"
#include "raytracingnew.cpp"
#include "redshift.cpp"
//#include "redshift_polish_doughnut.cpp" /*shafqat*/
//#include "polish_doughnut.cpp"          /*Shafqat*/
//#include "emission_angle_rth.cpp"       /*shafqat*/
#include "find_isco.cpp"                /*shafqat*/
#include "metric_rderivatives.cpp"      /*shafqat*/

#endif
