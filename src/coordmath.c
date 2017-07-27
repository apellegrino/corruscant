#include <math.h>

#define PI_OVER_180 (M_PI / 180.0)

/*
void radec2cart64(double * radec, double * xyz)
{
    double lon = radec[0] * PI_OVER_180;
    double lat = radec[1] * PI_OVER_180;
    xyz[0] = cos(lat) * cos(lon);
    xyz[1] = cos(lat) * sin(lon);
    xyz[2] = sin(lat);
}
*/

/* 
 * Projects RA, DEC points onto the unit sphere in Cartesian coordinates. Note
 * the transition from parallel arrays as input to an array of records on the
 * output.
 */
void radec2cart64(double * ra, double * dec, double * xyz, int n)
{
    int i;
    double lon, lat, coslat;

    for(i=0; i<n; i++) {
        lon = ra[i] * PI_OVER_180;
        lat = dec[i] * PI_OVER_180;
        coslat = cos(lat);
        xyz[3*i] = coslat * cos(lon);
        xyz[3*i+1] = coslat * sin(lon);
        xyz[3*i+2] = sin(lat);
    }
}

/* 
 * Converts spherical points (RA, DEC, DIST) to Cartesian coordinates. Note the
 * transition from parallel arrays as input to an array of records on the
 * output.
 */
void radecdist2cart64(double * ra, double * dec, double * dist, double * xyz,
                      int n)
{
    int i;
    double lon, lat, coslat;
    for(i=0; i<n; i++) {
        lon = ra[i] * PI_OVER_180;
        lat = dec[i] * PI_OVER_180;
        coslat = cos(lat);
        xyz[3*i] = dist[i] * coslat * cos(lon);
        xyz[3*i+1] = dist[i] * coslat * sin(lon);
        xyz[3*i+2] = dist[i] * sin(lat);
    }
}
