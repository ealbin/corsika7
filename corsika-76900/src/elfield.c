/*====================================================================*/

/*      elfield.c                                                     */

/*--------------------------------------------------------------------*/
/* HEADER DECLARATIONS & DEFINITIONS                                  */
#if defined(__GNUC__)&&__GNUC__ > 3
#include <stdlib.h>
#endif
#include <stdio.h>
#include <math.h>
/*--------------------------------------------------------------------*/
/* in C language: float: 32 bits,  range 3.4*1e-38  to 3.4*1e+38      */
/*                double: 64 bits, range 1.7*1e-308 to 1.7*1e+308     */
/*--------------------------------------------------------------------*/
/* DESCRIPTION of the procedure elfield.c:                            */
/*                                                                    */
/* INPUT:                                                             */
/*   double *x0: value of x-coordinate (cm)                           */
/*   double *y0: value of y-coordinate (cm)                           */
/*   double *z0: value of z-coordinate (cm)                           */
/*   NOTE: The x,y-cordinates are relative to the middle of the       */ 
/*                  lowest observation level.                         */
/*         The z-coordinate is relative to sea level                  */
/* OUTPUT:                                                            */
/*   double *ex: field strength (volt/cm) in x-direction              */
/*   double *ey: field strength (volt/cm) in y-direction              */
/*   double *ez: field strength (volt/cm) in z-direction              */
/*                                                                    */
/* For orientation of directions x, y, z see CORSIKA User's Guide     */
/*--------------------------------------------------------------------*/
/*-------- template for elfield.c subroutine -------------------------*/
/*--------------------------------------------------------------------*/
#ifdef __REDHAT__
void elfield__(double *x0, double *y0, double *z0,
               double *ex, double *ey, double *ez)
#else
void elfield_(double *x0, double *y0, double *z0,
              double *ex, double *ey, double *ez)
#endif
/*--------------------------------------------------------------------*/
{
/*----here the subroutine must be programmed -------------------------*/
 *ex=0.0;
 *ey=0.0;
 *ez=0.0;
}

