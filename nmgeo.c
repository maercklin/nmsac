/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMGEO.C - Functions for great-circle calculations
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, March 2006

Version: 2008-11-01

Note:
  All coordinates are in degrees and must be valid latitudes and longitudes.
  The signs of latitudes and longitudes are assumed to be positive in North 
  and East direction, and negative in South and West direction.

Credit:
  Some functions are based on the Aviation Formulary by Ed Williams
  (http://williams.best.vwh.net/avform.htm).
*************************************************************************/

/* Include files */
#include <stdlib.h>
#include <math.h>


/* Definitions */
#ifndef PI
#define PI 3.14159265358979322702
#endif
#define DDEGTORAD 0.01745329251994329547
#define DRADTODEG 57.2957795130823228646
#define WGS1984_A 6378137.0000
#define WGS1984_F 0.003352810665


/* Function prototypes */
double gc_azimuth(double lat1, double lon1, double lat2, double lon2);
double gc_delta(double lat1, double lon1, double lat2, double lon2);
void gc_inter(double *olat, double *olon, \
    double lat1, double lon1, double lat2, double lon2, double frac);
double gc_dist(double lat1, double lon1, double lat2, double lon2,
    double a, double f);



/************************************************************************/
/* Functions for great-circle calculations                              */
/************************************************************************/

double gc_azimuth(double lat1, double lon1, double lat2, double lon2) 
/************************************************************************
gc_azimuth - compute the azimuth of two points on a sphere

Input:
lat1    latitude  of point1 in degrees
lon1    longitude of point1 in degrees
lat2    latitude  of point2 in degrees
lon2    longitude of point2 in degrees

Output: returns azimuth in degrees

Notes:
Azimuths are measured from North in clockwise direction 0-359.999 degrees.
Reverse point1 and point2 to compute the backazimuth.
The signs of input latitudes and longitudes are positive in North and East
direction and negative in South and West direction.

Author: Nils Maercklin, 2006-03-19
*************************************************************************/
{
    /* Output azimuth value */
    double azimuth;

    /* Convert input from degrees to radians */
    lat1 *= DDEGTORAD;
    lon1 *= DDEGTORAD;
    lat2 *= DDEGTORAD;
    lon2 *= DDEGTORAD;

    /* Compute azimuth */
    if (cos(lat1)<0.000000001) {
        azimuth = (lat1>0.0) ? PI : 0.0;  /* pole */
    }
    else {
        azimuth = atan2(sin(lon2-lon1)*cos(lat2), \
            cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1));
        azimuth = fmod(2.0*PI + azimuth, 2.0*PI);
    }

    /* Return azimuth in degrees */
    return azimuth * DRADTODEG;
}



double gc_delta(double lat1, double lon1, double lat2, double lon2)
/************************************************************************
gc_delta - distance in degrees between two points on a sphere

Input:
lat1    latitude  of point1 in degrees
lon1    longitude of point1 in degrees
lat2    latitude  of point2 in degrees
lon2    longitude of point2 in degrees

Output: returns distance in degrees

Author: Nils Maercklin, 2006-03-19
*************************************************************************/
{
    /* Computed distance in degrees */
    double d;

    /* Convert input values from degrees to radians */
    lat1*=DDEGTORAD; lon1*=DDEGTORAD;
    lat2*=DDEGTORAD; lon2*=DDEGTORAD;

    /* Compute the distance on the sphere */
    d = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));

    /* Return distance in degrees */
    return d*DRADTODEG;
}



void gc_inter(double *olat, double *olon, \
    double lat1, double lon1, double lat2, double lon2, double frac)
/************************************************************************
gc_inter - intermediate point on a great-circle path 

Input:
lat1    latitude  of point1 in degrees
lon1    longitude of point1 in degrees
lat2    latitude  of point2 in degrees
lon2    longitude of point2 in degrees
frac    fraction of distance for output intermediate point [0...1]

Output:
olat    latitude of intermediate point in degrees
olat    longitude of intermediate point in degrees

Note:
The two input points must not be antipodes (great-circle path undefined).
An input frac=0 yields point1 and frac=1 point2, respectively.

Author: Nils Maercklin, 2006-03-19
*************************************************************************/
{
    /* Internal variables */
    double a, b, d, x, y, z;

    /* Convert input from degrees to rad */
    lat1*=DDEGTORAD; lon1*=DDEGTORAD;
    lat2*=DDEGTORAD; lon2*=DDEGTORAD;

    /* Distance between points */
    d = DDEGTORAD * gc_delta(lat1, lon1, lat2, lon2);

    /* Some convenient abbreviations */
    a = sin((1.0-frac)*d)/sin(d);
    b = sin(frac*d)/sin(d);

    x = a*cos(lat1)*cos(lon1) +  b*cos(lat2)*cos(lon2);
    y = a*cos(lat1)*sin(lon1) +  b*cos(lat2)*sin(lon2);
    z = a*sin(lat1)           +  b*sin(lat2);

    /* Output computed coordinates in degrees */
    *olat = atan2(z,sqrt(x*x + y*y)) * DRADTODEG;
    *olon = atan2(y,x) * DRADTODEG;
}



double gc_dist(double lat1, double lon1, double lat2, double lon2,
    double a, double f)
/************************************************************************
gc_dist - distance between two points along great circle path 

Input:
lat1    latitude  of point1 in degrees
lon1    longitude of point1 in degrees
lat2    latitude  of point2 in degrees
lon2    longitude of point2 in degrees
a       semi-major axis of ellipsoid in units of length
f       flattening of ellipsoid (0 for a sphere)

Output: returns distance in units of length

Notes:
The distance on an ellipsoid is computed via an approximate algorithm, 
translated from a GAWK function (NM, 2004; adopted from a C++ method 
published at www.codeguru.com). 

Author: Nils Maercklin, 2006-03-19
*************************************************************************/
{
    /* Internal Variables */
    double edist;
    double dd,ff,gg,ll,ss,cc,ww,rr,hh1,hh2,ss1,cc2;

    /* Approximate distance on ellipsoid */
    if (f) {
        /* Convert coordinates from degrees to radians */
        lat1*=DDEGTORAD; lon1*=DDEGTORAD;
        lat2*=DDEGTORAD; lon2*=DDEGTORAD;

        /* Definition of several abbreviations */
        ff = (lat1 + lat2) / 2.0;
        gg = (lat1 - lat2) / 2.0;
        ll = (lon1 - lon2) / 2.0;

        ss = sin(gg)*sin(gg)*cos(ll)*cos(ll) + cos(ff)*cos(ff)*sin(ll)*sin(ll);
        cc = cos(gg)*cos(gg)*cos(ll)*cos(ll) + sin(ff)*sin(ff)*sin(ll)*sin(ll);
        ww = atan2(sqrt(ss),sqrt(cc));

        if (ww==0.) { return 0.0; } /* quick hack (NM, 2004-10-01)*/
        rr = sqrt(ss*cc)/ww;

        hh1 = (3.0*rr - 1.0) / (2.0*cc);
        hh2 = (3.0*rr + 1.0) / (2.0*ss);
        dd  = 2.0*ww*a;


        ss1 = sin(ff)*sin(ff)*cos(gg)*cos(gg);
        cc2 = cos(ff)*cos(ff)*sin(gg)*sin(gg);

        /* Compute distance on ellipsoid */
        edist = dd * (1.0 + f*hh1*ss1 - f*hh2*cc2);
    }
    /* Flattening f=0, i.e. "exact" computation on a sphere */
    else {
        edist = a * DDEGTORAD * gc_delta(lat1,lon1,lat2,lon2);
    }

    /* Return distance in units of length */
    return edist;
}

/* END OF FILE */
