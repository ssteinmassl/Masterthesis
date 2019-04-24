/***************************************
*
*                     FILE  STAR1.C
*
*   Functions concerned with the compact star ("star 1" or the
*   "primary star") are in this file.
*
****************************************/

#include "header.h"

double Star1TotFlux( double distance )
/*****************************************************
*
*   This function returns the integrated flux from star 1.
*   The flux has been integrated over wavelength and over 
*   the surface of the star.
*
**************************************************************/
{
  double totalflux;

  totalflux = star1.L / (FOURPI * distance * distance );
  return( totalflux );
}


double Star1Flambda( char *filter, double minlambda, double maxlambda )
/*****************************************************
*
*   This function returns the contribution of star 1 to the
*   observed spectrum at wavelength lambda.  Note that the wavelengths
*   must be in centimeters but the returned flux is in
*   ergs/sec/cm^2/Angstrom.
*
**************************************************************/
{
   double flux;

   if( strcmp( filter, "SQUARE") == 0 ) {
      flux = ( star1.L / (4.0 * star1.sigmaT4) )
              * BBSquareIntensity( star1.T, minlambda, maxlambda );
   }
   else {
      flux = ( star1.L / (4.0 * star1.sigmaT4) )
	      * BBFilterIntensity( star1.T, filter );
   }

   return( flux );
}
