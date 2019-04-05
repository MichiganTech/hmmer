/* drand48.c
 * 
 * Portable random number generator, and sampling routines.
 *
 * SRE, Tue Oct  1 15:24:11 2002 [St. Louis]
 * CVS $Id: drand48.c,v 1.1 2002/10/09 14:26:09 eddy Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sre_random.h"


/* Function: sre_srandom()
 * 
 * Purpose:  Initialize with a random seed. Seed must be
 *           >= 0 to work; we silently enforce this.
 */
// void
// sre_srandom(int seed)
// {
//   if (seed < 0)  seed = -1 * seed;
//   if (seed == 0) seed = 42;
//   sre_randseed = seed;
// }


/* Function: drand48_positive()
 *
 * Purpose:  Assure 0 < x < 1 (positive uniform deviate)
 */
double
drand48_positive(
){
  double x;
  do { x = drand48(); } while (x == 0.0);
  return x;
}

/* Function: ExponentialRandom()
 *
 * Purpose:  Pick an exponentially distributed random variable
 *           0 > x >= infinity
 *           
 * Args:     (void)
 *
 * Returns:  x
 */
double
ExponentialRandom(
){
  double x;

  do x = drand48(); while (x == 0.0);
  return -log(x);
}    

/* Function: Gaussrandom()
 * 
 * Pick a Gaussian-distributed random variable
 * with some mean and standard deviation, and
 * return it.
 * 
 * Based on RANLIB.c public domain implementation.
 * Thanks to the authors, Barry W. Brown and James Lovato,
 * University of Texas, M.D. Anderson Cancer Center, Houston TX.
 * Their implementation is from Ahrens and Dieter, "Extensions 
 * of Forsythe's method for random sampling from the normal
 * distribution", Math. Comput. 27:927-937 (1973).
 *
 * Impenetrability of the code is to be blamed on its FORTRAN/f2c lineage.
 * 
 */
double
Gaussrandom(double mean, double stddev)
{
  static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
  };
  static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
  };
  static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
  };
  static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
  };
  static long i;
  static double snorm,u,s,ustar,aa,w,y,tt;

  u = drand48();
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
   * START CENTER
   */
  ustar = u-(double)i;
  aa = *(a+i-1);
S40:
  if(ustar <= *(t+i-1)) goto S60;
  w = (ustar-*(t+i-1))**(h+i-1);
S50:
  /*
   * EXIT   (BOTH CASES)
   */
  y = aa+w;
  snorm = y;
  if(s == 1.0) snorm = -y;
  return (stddev*snorm + mean);
S60:
  /*
   * CENTER CONTINUED
   */
  u = drand48();
  w = u*(*(a+i)-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
S70:
  tt = u;
  ustar = drand48();
S80:
  if(ustar > tt) goto S50;
  u = drand48();
  if(ustar >= u) goto S70;
  ustar = drand48();
  goto S40;
S100:
  /*
   * START TAIL
   */
  i = 6;
  aa = *(a+31);
  goto S120;
S110:
  aa += *(d+i-1);
  i += 1;
S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
S140:
  w = u**(d+i-1);
  tt = (0.5*w+aa)*w;
  goto S160;
S150:
  tt = u;
S160:
  ustar = drand48();
  if(ustar > tt) goto S50;
  u = drand48();
  if(ustar >= u) goto S150;
  u = drand48();
  goto S140;
}
