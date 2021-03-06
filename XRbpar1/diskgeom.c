/***************************************
*
*                     FILE  DISKGEOM.C
*
*   Functions concerned with the geometry of the accretion disk 
*   around the compact star are in this file.
*
****************************************/

#include "header.h"

void MakeDiskTiles( void )
/************************************************
*
*   This function distributes tiles over the surface of the disk
*   and calculates various properties of the tiles.  The tiles
*   all have roughly the same area.
*
*   The first tile is tile #1 (not zero)
*
**************************************************/
{
   FILE *out;
   char filename[40];

   long i, nringMain, nringEdge, ring, nzeta, tilenumber,
        maintiles, edgetiles, band, mintiles;
   double maindiskarea, rimheight, edgearea, torusheight, extra,
          extraAtorus, extraRtorus, extraArim, extraRrim, flatArea, 
          rougharea, targetarea, targetda, daMain, daRim;
   double x, a, a1, a2, h, h1, h2, slope, projection, sum1, sum2, 
          meanh1, meanh2, ringarea, zeta, zeta1, zeta2, rho1, rho2,
          dzeta, height, dh;

   if( strcmp( verbose, "ON") == 0 ) {
     printf(" Begin making disk tiles.\n");
   }

   /************************************************************
   *
   *   Calculate the typical area and dimensions of the tiles.
   *
   *   Method:  Calculate a rough area for the disk and divide it
   *       by the target number of tiles to get a target area for
   *       the tiles.  The target edge length of the tiles is 
   *           targetda  =  sqrt( target area).
   *
   **************************************************************/
   maindiskarea = PI * ( maindisk.amax * maindisk.amax
                              - maindisk.amin * maindisk.amin);
   rimheight= 0.0;
   if( strcmp( control.diskrim, "ON") == 0) {
      rimheight = 0.5 * (diskrim.Hmin + diskrim.Hmax);
   }
   edgearea = 2.0 * PI * maindisk.amax * (maindisk.Hmax + rimheight);
   extraArim = 0.0;
   extraRrim = 0.0;
   if( strcmp( control.diskrim, "ON" ) == 0 ) {
      flatArea = PI * maindisk.amax * maindisk.amax
                         - PI * (maindisk.amax - diskrim.awidth) 
                              * (maindisk.amax - diskrim.awidth);
      extra = 0.5 * (rimheight * rimheight)
                                 / (diskrim.awidth * diskrim.awidth);
      extraArim = flatArea * extra;
      extraRrim = diskrim.awidth * extra;
   }
   extraAtorus = 0.0;
   extraRtorus = 0.0;
   if( strcmp( control.disktorus, "ON" ) == 0 ) {
      flatArea =   PI * (disktorus.azero + 0.5 * disktorus.awidth)
                      * (disktorus.azero + 0.5 * disktorus.awidth)
                 - PI * (disktorus.azero - 0.5 * disktorus.awidth)
                      * (disktorus.azero - 0.5 * disktorus.awidth);
      torusheight = 0.5 * (disktorus.Hmax + disktorus.Hmin);
      extra =  0.5 * (torusheight * torusheight)
                      / (0.5 * disktorus.awidth * 0.5 * disktorus.awidth);
      extraAtorus = flatArea * extra;
      extraRtorus = disktorus.awidth * extra;
   }
   rougharea = 2.0 * ( maindiskarea +  edgearea + extraAtorus + extraArim);

   targetarea = rougharea / disk.targetNtiles;
   targetda = sqrt( targetarea );
   nringMain = 0.5 + (maindisk.amax - maindisk.amin
                              + extraRtorus + extraRrim) / targetda;
   if( nringMain < 2 )
      Quit("disk.targetNtiles too small to cover the disk properly.");
   daMain = (maindisk.amax - maindisk.amin - 1.0e-7*maindisk.amax) / nringMain;
   nringEdge = 0.5 + (maindisk.Hmax + rimheight) / targetda;
   if( nringEdge < 1 )
      nringEdge = 1;

   /******************************************************
   *
   *   Now lay out the tiles.  Note that the tiles are first
   *   distributed in (a,zeta) space and then converted to
   *   (rho,zeta) space.
   *   Every ring must have at least mintiles tiles per ring.
   *
   ***********************************************************/

   /*********************************************************
   *
   *   Tiles on the top surface of the disk.
   *
   *********************************************************/

   mintiles = 30;
   tilenumber = 0;
   for( ring = 1; ring <= nringMain; ring++) {
      a1 = (1.0 + 1.0e-7) * maindisk.amin + (ring -1) * daMain;
      a2 = a1 + daMain;
      h1 = MainDiskH( a1 );
      h2 = MainDiskH( a2 );
      a = 0.5 * (a1 + a2);
      if( strcmp( control.disktorus, "ON") == 0 ) {
	 if( (a > (disktorus.azero - 0.5 * disktorus.awidth))
	              && (a < (disktorus.azero + 0.5 * disktorus.awidth)) ) {
            h1 = DiskTopH( a1, disktorus.ZetaHmax );
            h2 = DiskTopH( a2, disktorus.ZetaHmax );
	 }
      }
      if( strcmp( control.diskrim, "ON") == 0 ) {
	 if( a > (maindisk.amax - diskrim.awidth) ) {
            h1 = DiskTopH( a1, diskrim.ZetaHmax );
            h2 = DiskTopH( a2, diskrim.ZetaHmax );
	 }
      }
      slope = (h2 - h1) / (a2 - a1);
      projection = sqrt( 1.0 + slope*slope );
      if( projection > 3.0 ) 
	 projection = 4.0;
      ringarea = projection * ( PI*a2*a2 - PI*a1*a1);
      nzeta = 0.5 + (ringarea / targetarea);
      if( nzeta < mintiles )
         nzeta = mintiles;
      dzeta = TWOPI / nzeta;
      for( i = 1; i <= nzeta; i++) {
         tilenumber += 1;
         if( tilenumber >= (MAXDISKTILES / 2) )
	    Quit("Attempted to make too many disk tiles.");
         zeta1 = (i - 1) * dzeta;
         zeta2 = zeta1 + dzeta - 1.0e-10;
         zeta = 0.5 * (zeta1 + zeta2);
         TDiska[tilenumber] = a;
         TDiskZeta[tilenumber] = zeta;
         TDiskH[tilenumber] = 0.5 * ( DiskTopH( a1, zeta) 
                                   +  DiskTopH( a2, zeta) );
         TDiskRho[tilenumber] = AToRho( a, zeta );
         TDiskx[tilenumber] = TDiskRho[tilenumber] * sin( zeta );
	 TDisky[tilenumber] = TDiskH[tilenumber];
         TDiskz[tilenumber] = TDiskRho[tilenumber] * cos( zeta ) + syspars.a;
         TDisknormCyl[tilenumber] = DiskNormal( "TOP", a1, a2, zeta1, zeta2);
         TDisknormCart[tilenumber] = Cyl2Cart( TDisknormCyl[tilenumber], 
                                          zeta );
         TDiskdS[tilenumber] = TopTileArea( tilenumber, a1, a2, zeta1, zeta2);
         TDiskT[tilenumber] = DiskTopT( a, zeta );
         TDiskT4[tilenumber] = pow( TDiskT[tilenumber], 4.0);
      }
   }
   maintiles = tilenumber;

   /**************************************************
   *
   *   Now put the tiles on the edge of the disk.  The tiles
   *   are in a different order than on the surface.
   *
   ***************************************************/
   nzeta = 0.5 + edgearea / (nringEdge * targetarea);
   if( nzeta < mintiles )
      nzeta = mintiles;
   dzeta = TWOPI / nzeta;
   for( i = 1; i <= nzeta; i++) {
      zeta1 = (i - 1) * dzeta;
      zeta2 = zeta1 + dzeta - 1.0e-10;
      zeta = 0.5 * (zeta1 + zeta2);
      height = maindisk.Hmax + DiskRimH( maindisk.amax, zeta );
      dh = height / nringEdge;
      for( ring = 1; ring <= nringEdge; ring++) {
         tilenumber += 1;
         if( tilenumber >= (MAXDISKTILES / 2) )
	    Quit("Attempted to make too many disk tiles.");
         TDiska[tilenumber] = maindisk.amax;
         h = (ring - 0.5) * dh;
         TDiskRho[tilenumber] = AToRho( maindisk.amax, zeta );
         TDiskZeta[tilenumber] = zeta;
         TDiskH[tilenumber] = h;
         TDiskx[tilenumber] = TDiskRho[tilenumber] * sin( zeta );
	 TDisky[tilenumber] = TDiskH[tilenumber];
         TDiskz[tilenumber] = TDiskRho[tilenumber] * cos( zeta ) + syspars.a;
         TDisknormCyl[tilenumber] = DiskNormal( "EDGE", a1, a2, zeta1, zeta2);
         TDisknormCart[tilenumber] = Cyl2Cart( TDisknormCyl[tilenumber],
                                                      zeta );
         TDiskdS[tilenumber] = EdgeTileArea( tilenumber, zeta1, zeta2, dh);
         TDiskT[tilenumber] = DiskEdgeT( zeta );
         TDiskT4[tilenumber] = pow( TDiskT[tilenumber], 4.0);
      }
   }
   edgetiles = tilenumber - maintiles;
   
   /************************************************************
   *
   *   Duplicate the tiles on the bottom of the disk.  This is 
   *   necessary because the bottom tiles can irradiate and heat 
   *   parts of the secondary star that are visible.
   *
   *************************************************************/
   for( i = 1; i <= tilenumber; i++) {
      if( (i+tilenumber) >= (MAXDISKTILES / 2) )
	    Quit("Attempted to make too many disk tiles.");
      TDiska[i+tilenumber]    =  TDiska[i];
      TDiskRho[i+tilenumber]  =  TDiskRho[i];
      TDiskZeta[i+tilenumber] =  TDiskZeta[i];
      TDiskH[i+tilenumber]    = -TDiskH[i];
      TDiskx[i+tilenumber]    =  TDiskx[i];
      TDisky[i+tilenumber]    = -TDisky[i];
      TDiskz[i+tilenumber]    =  TDiskz[i];
      TDiskdS[i+tilenumber]   =  TDiskdS[i];
      TDisknormCyl[i+tilenumber].rho  =  TDisknormCyl[i].rho;
      TDisknormCyl[i+tilenumber].zeta =  TDisknormCyl[i].zeta;
      TDisknormCyl[i+tilenumber].h    = -TDisknormCyl[i].h;
      TDisknormCart[i+tilenumber].x   =  TDisknormCart[i].x;
      TDisknormCart[i+tilenumber].y   = -TDisknormCart[i].y;
      TDisknormCart[i+tilenumber].z   =  TDisknormCart[i].z;
      TDiskT[i+tilenumber]  = TDiskT[i];
      TDiskT4[i+tilenumber] = TDiskT4[i];
   }
   disk.Ntiles = 2 * tilenumber;

   /************************************************************
   *
   *   The last step:  Calculate the specific intensity emitted
   *   by each disk tile.  The disk is assumed to emit angle-independent
   *   black body intensities.  The mean intensities depend, then,
   *   on the bandpass:
   *
   *************************************************************/
   for( band = 1; band <= orbit.nbands; band++) {
      if( strcmp( orbit.filter[band], "SQUARE") == 0 ) {
         for( i = 1; i <= disk.Ntiles; i++) {
            TDiskI[band][i] = BBSquareIntensity( TDiskT[i], 
                                  orbit.minlambda[band],
                                  orbit.maxlambda[band]);
         }
      }
      else {
         for( i = 1; i <= disk.Ntiles; i++) {
            TDiskI[band][i] = BBFilterIntensity( 
                                  TDiskT[i], orbit.filter[band]);
         }
      }
   }

   if( strcmp( control.diagnostics, "INSPECTDISKTILES") == 0 ) {
      InspectDiskTiles( targetarea, nringMain, maintiles, 
                            nringEdge, edgetiles);
      Quit("Quit after INSPECTDISKTILES.");
   }

   return;
}


double DiskTopH( double a, double zeta )
/***********************************************
*
*   This function returns the upper boundary of the disk
*   at position (rho, zeta).
*
*******************************************************/
{
   double hdisk, hrim, htorus, height;

   hdisk = MainDiskH( a );
   hrim  = DiskRimH( a, zeta );
   htorus = DiskTorusH( a, zeta );
   height = hdisk + hrim + htorus;

   return( height );
}


double DiskBottomH( double a, double zeta )
/***********************************************
*
*   This function returns the lower boundary of the disk
*   at position (rho, zeta).
*   The current version assumes that the disk is symmetric
*   about the orbital plane.
*
*******************************************************/
{
   double hdisk, hrim, htorus, height;

   hdisk = MainDiskH( a );
   hrim  = DiskRimH( a, zeta );
   htorus = DiskTorusH( a, zeta );
   height = -( hdisk + hrim + htorus );

   return( height );
}


double MainDiskH( double a )
/**************************************
*
*   This function returns the height of the main, axi-symmetric
*   part of the disk.  The height has the functional form
*       h = Hmax * ( (a - amin) / (amax - amin) )^Hpow
*   for
*       amin < a < amax
*   and
*       h = 0 otherwise.
*   Note:  To make a disk with constant thickness equal to
*          Hmax, just set Hpow = 0.
*
*   The weird factor (1.0 + 1.0e-7) is protection against 
*   roundoff error at the edge of the disk.
*
*****************************************/
{
   double x, height;

   if( a < maindisk.amin ) {
      height = 0.0;
   }
   else if ( a > ( (1.0 + 1.0e-7) * maindisk.amax ) ) {
      height = 0.0;
   }
   else {
      x = (a - maindisk.amin) / (maindisk.amax - maindisk.amin);
      if( x == 0.0 ) {
	 height = maindisk.Hmax;
      }
      else {
         height = maindisk.Hmax * pow( x, maindisk.Hpow );
      }
   }
   return( height );
}


double DiskRimH( double a, double zeta )
/**************************************
*
*   This function returns the height of the disk rim at 
*   position (a, zeta) .  The cross section of the rim is 
*   modeled as a truncated ellipse added on top of the main 
*   disk with
*      vertical semi-major axis   = H
*      horizontal semi-minor axis = diskrim.awidth
*
*       (h/H)^2 + ( (maindisk.amax -a) / awidth )^2 = 1
*
*   The rim height H is a function of zeta but its width is not.
*   There are currently two possibilities for the zeta dependence
*
*   1) SINUSOID
*      H = 0.5 * (Hmax + Hmin)
*            + 0.5* (Hmax - Hmin) * cos( zeta - zetaHmax );
*
*   2) POINT
*      The disk rim height and temperature is defined by a set 
*      of points, one point per DISKRIMPARS= line in the parameter file:
*         DISKRIMPARS=  POINT   Zeta1   H1   T1
*         DISKRIMPARS=  POINT   Zeta2   H2   T2    
*            .        .       .     .    .
*            .        .       .     .    .
*      The Zetas need not be in order.
*      The heights are linearly interpolated between the specified points.
*
*   The weird factor (1.0 + 1.0e-7) is protection against 
*   roundoff error at the edge of the disk.
*
*****************************************/
{
   long i;
   double Hzeta, zetalow, zetahigh, Hlow, Hhigh, slope, 
          x, y, height;

   if( zeta > TWOPI )
      Quit("zeta greater than TWOPI in DiskRimH.");
   if( zeta < 0.0 )
      Quit("zeta less than zero in DiskRimH.");

   if( a < (maindisk.amax - diskrim.awidth) ) {
      height = 0.0;
      return( height );
   } 
   if( a > ( (1.0 + 1.0e-7) * maindisk.amax ) ) {
      height = 0.0;
      return( height );
   }

   if( strcmp( diskrim.type, "SINUSOID" ) == 0 ) {
      Hzeta =    0.5 * ( diskrim.Hmax + diskrim.Hmin )
               + 0.5 * ( diskrim.Hmax - diskrim.Hmin ) 
                           * cos( zeta - diskrim.ZetaHmax );
   }
   else if( strcmp( diskrim.type, "POINT" ) == 0 ) {
      if( diskrim.points == 1 ) {
         Hzeta = diskrim.PointH[1];
      }
      else {
	 for( i = 1; i <= diskrim.points; i++) {
	    zetalow = diskrim.PointZeta[i];
            Hlow = diskrim.PointH[i];
            if( i < diskrim.points ) {
	       zetahigh = diskrim.PointZeta[i+1];
               Hhigh = diskrim.PointH[i+1];
            }
            else {
	       zetahigh = TWOPI;
               Hhigh = diskrim.PointH[1];
            }
	    if( (zeta >= zetalow) && (zeta < zetahigh) ) {
	       slope = (Hhigh - Hlow) / (zetahigh - zetalow);
               Hzeta = Hlow + slope * (zeta - zetalow);
               break;
            }
	 }
      }
   }
   else
      Quit("Unrecognized disk rim type in DiskRimH.");

   x = (maindisk.amax - a) / diskrim.awidth;
   y = sqrt( 1.0 - x*x);
   height = Hzeta * y;

   return( height );
}


double DiskTorusH( double a, double zeta )
/**************************************
*
*   This function returns the height of the disk torus
*   at position (a, zeta).   The cross section of the torus is 
*   modeled as an ellipse added on top of the main disk with
*      ellipse center             = disktorus.azero
*      vertical semi-major axis   = H
*      horizontal semi-minor axis = 0.5 * disktorus.awidth
*
*       (h/H)^2 + ( (disktorus.azero -a) / awidth )^2 = 1
*
*   The torus height H is a function of zeta but its width is not. 
*   There are currently two possibilities for the zeta dependence
*
*   1) SINUSOID
*      H = 0.5 * (Hmax + Hmin)
*            + 0.5* (Hmax - Hmin) * cos( zeta - zetaHmax );
*
*   2) POINT
*      The disk torus height and temperature is defined by a set of
*      points, one point per DISKTORUSPARS= line in the parameter file:
*         DISKTORUSPARS=  POINT   Zeta1   H1   T1
*         DISKTORUSPARS=  POINT   Zeta2   H2   T2    
*            .        .       .     .    .
*            .        .       .     .    .
*      The heights are linearly interpolated between the specified points.
*
*      The Zetas must be in increasing order and disktorus.PointZeta[1]
*      must be 0 degrees (this avoids messy computer code).
*      The heights are linearly interpolated between the specified points.
*
*****************************************/
{
   long i;
   double Hzeta, zetalow, zetahigh, Hlow, Hhigh, slope, 
          x, y, height;

   if( zeta > TWOPI )
      Quit("zeta greater than TWOPI in DiskTorusH.");
   if( zeta < 0.0 )
      Quit("zeta less than zero in DiskTorusH.");

   if( a < (disktorus.azero - 0.5 * disktorus.awidth) ) {
      height = 0.0;
      return( height );
   } 
   if( a > (disktorus.azero + 0.5 * disktorus.awidth) ) {
      height = 0.0;
      return( height );
   }

   if( strcmp( disktorus.type, "SINUSOID" ) == 0 ) {
      Hzeta =    0.5 * ( disktorus.Hmax + disktorus.Hmin )
               + 0.5 * ( disktorus.Hmax - disktorus.Hmin ) 
                           * cos( zeta - disktorus.ZetaHmax );
   }
   else if( strcmp( disktorus.type, "POINT" ) == 0 ) {
      if( disktorus.points == 1 ) {
         Hzeta = disktorus.PointH[1];
      }
      else {
	 for( i = 1; i <= disktorus.points; i++) {
	    zetalow = disktorus.PointZeta[i];
            Hlow = disktorus.PointH[i];
            if( i < disktorus.points ) {
	       zetahigh = disktorus.PointZeta[i+1];
               Hhigh = disktorus.PointH[i+1];
            }
            else {
	       zetahigh = TWOPI;
               Hhigh = disktorus.PointH[1];
            }
	    if( (zeta >= zetalow) && (zeta < zetahigh) ) {
	       slope = (Hhigh - Hlow) / (zetahigh - zetalow);
               Hzeta = Hlow + slope * (zeta - zetalow);
               break;
            }
	 }
      }
   }
   else
      Quit("Unrecognized disk torus type in DiskTorusH.");

   x = (disktorus.azero - a) / (0.5 * disktorus.awidth);
   y = sqrt( 1.0 - x*x);
   height = Hzeta * y;

   return( height );
}


struct CylVector DiskNormal( char *where, double a1, double a2,
                                  double zeta1, double zeta2)
/***************************************************
*
*   This function calculates the unit vector that is
*   normal to the surface of a disk tile.  Note: the vector is
*   normal to the tile, not normal to the surface of the disk
*   underlying the tile.
*
*****************************************************/
{
   struct CylVector norm;
   double a, zeta, rho1, rho2, rho, h1, h2, DhDrho,
          DhDzeta, DrhoDzeta, beta;

   if( (a2 - a1) <= 0.0 )
      Quit("a2 le a1 in DiskNormal().");
   if( (zeta2 - zeta1) <= 0.0 )
      Quit("zeta2 le zeta1 in DiskNormal().");
   a = 0.5 * (a1 + a2);
   zeta = 0.5 * (zeta1 + zeta2);
   if( strcmp( where, "TOP") == 0 ) {
      rho1 = AToRho( a1, zeta);
      rho2 = AToRho( a2, zeta);
      rho =  AToRho( a,  zeta);
      h1 = DiskTopH( a1, zeta);
      h2 = DiskTopH( a2, zeta);
      DhDrho = (h2 - h1) / (rho2 - rho1);
      h1 = DiskTopH( a, zeta1);
      h2 = DiskTopH( a, zeta2);
      DhDzeta = (h2 - h1) / (zeta2 - zeta1);
      norm.rho  = -DhDrho;
      norm.zeta = -(1.0 / rho) * DhDzeta;
      norm.h    =  1.0;
    }
   else if( strcmp( where, "EDGE") == 0 ) {
      rho1 = AToRho( a, zeta1);
      rho2 = AToRho( a, zeta2);
      rho =  AToRho( a, zeta);
      DrhoDzeta = (rho2 - rho1) / (zeta2 - zeta1);
      norm.rho  = 1.0;
      norm.zeta = -(1.0 / rho) * DrhoDzeta;
      norm.h    = 0.0;
   }
   else {
      Quit("Unrecognized location in DiskNormal.");
   }

   beta = 1.0 / sqrt( norm.rho*norm.rho 
                         + norm.zeta*norm.zeta + norm.h*norm.h);
   norm.rho  = beta * norm.rho;
   norm.zeta = beta * norm.zeta;
   norm.h    = beta * norm.h;

   return( norm );
}


double TopTileArea( long itile, double a1, double a2, 
                                double zeta1, double zeta2)
/**********************************************************
*
*   This function calculates the area of tiles on the top surface
*   of the disk extending from a1 to a2, and from zeta1 to zeta2.
*   Note that the tile is finite in size, so the expression
*   for flatarea is not just  d S = rho d rho d zeta.
*
*****************************************************************/
{
   double zeta, dzeta, rho1, rho2, flatarea, Area;

   if( zeta1 < 0.0 )
      Quit("zeta less than zero in TopTileArea.");
   if( zeta2 > TWOPI )
      Quit("zeta greater than TWOPI in TopTileArea.");
   zeta = 0.5 * (zeta1 + zeta2);
   dzeta = zeta2 - zeta1;
   if( dzeta < 0.0 )
      Quit("dzeta is less than zero in TopTileArea.");

   rho1 = AToRho( a1, zeta);
   rho2 = AToRho( a2, zeta);
   flatarea = PI * (rho2*rho2 - rho1*rho1) * ( dzeta / TWOPI );;

   if( TDisknormCyl[itile].h == 0.0 )
      Quit("TDisknormCyl.h equals zero in TopTileArea.");
   Area = flatarea / TDisknormCyl[itile].h;

   return( Area );
}


double EdgeTileArea( long itile, double zeta1, double zeta2, double dh)
/**********************************************************
*
*   This function calculates the area of edge tiles.
*
*****************************************************************/
{
   double zeta, dzeta, flatarea, Area;

   if( zeta1 < 0.0 )
      Quit("zeta less than zero in EdgeTileArea.");
   if( zeta2 > TWOPI )
      Quit("zeta greater than TWOPI in EdgeTileArea.");
   zeta = 0.5 * (zeta1 + zeta2);
   dzeta = zeta2 - zeta1;
   if( dzeta < 0.0 )
      Quit("dzeta less than zero in EdgeTileArea.");

   flatarea = TDiskRho[itile] * dzeta * dh;

   if( TDisknormCyl[itile].rho == 0.0 )
      Quit("TDisknormCyl.rho equals zero in EdgeTileArea.");
   Area = flatarea / TDisknormCyl[itile].rho;

   return( Area );
}


double RhoToA( double rho, double zeta )
/************************************************
*
*   Calculates a from (rho, zeta).
*
**************************************************/
{
   double a;

   if( disk.e == 0.0 ) {
      a = rho;
   }
   else {
      a = rho * (1.0 + disk.e * cos(zeta - disk.zetazero) ) 
                              / (1.0 - disk.e * disk.e);
   }

   return( a );
}


double AToRho( double a, double zeta )
/************************************************
*
*   Calculates rho from ( a, zeta).
*
**************************************************/
{
   double rho;

   if( disk.e == 0.0 ) {
      rho = a;
   }
   else {
      rho =  ( a * (1.0 - disk.e * disk.e) ) 
                   / (1.0 + disk.e * cos(zeta - disk.zetazero) );
   }
                             
   return( rho );
}


void MakeDiskRim( void )
/******************************************************
*
*   This function is invoked if the disk rim is ON and the rim
*   is specified by a set of discrete points along its rim.  
*   It does the following:
*      1)  Sort the rim points in order of increasing zeta.
*      2)  Check for duplicate points.
*      3)  Create a point a zeta = 0 if one does not already
*            exist by linearly interpolating in zeta.
*      4)  Finds maximum and minimum values for H and T.
*
*   Note:  The function assumes that zeta is given in degrees.
*
*************************************************************/
{
   long i, j;
   double dummy, zeta2, zetaN;

   if( diskrim.points <= 0 )
      Quit("SortRimPoints:  diskrim.points le 0.");
   else if( diskrim.points > 1 ) {
      for( i = 1; i < diskrim.points; i++ ) {
         for( j = i+1; j <= diskrim.points; j++) {
            if( diskrim.PointZeta[i] > diskrim.PointZeta[j] ) {
 	       dummy = diskrim.PointZeta[i];
               diskrim.PointZeta[i] = diskrim.PointZeta[j];
               diskrim.PointZeta[j] = dummy;
               dummy = diskrim.PointH[i];
               diskrim.PointH[i] = diskrim.PointH[j];
               diskrim.PointH[j] = dummy;
               dummy = diskrim.PointT[i];
               diskrim.PointT[i] = diskrim.PointT[j];
               diskrim.PointT[j] = dummy;
            }
         }
      }
   }

   if( diskrim.points > 1 ) {
      for( i = 1; i <= diskrim.points; i++) {
         if( diskrim.PointZeta[i] == diskrim.PointZeta[i+1] )
            Quit("DISKRIMPARS: At least two points have the same PointZeta.");
      }
   }
   if( diskrim.PointZeta[1] != 0.0 ) {
      diskrim.points += 1;
      for( i = diskrim.points; i >= 2; i--) {
         diskrim.PointZeta[i] = diskrim.PointZeta[i-1];
         diskrim.PointH[i]    = diskrim.PointH[i-1];
         diskrim.PointT[i]    = diskrim.PointT[i-1];
      }
      diskrim.PointZeta[1] = 0.0;
      zeta2 = diskrim.PointZeta[2];
      zetaN = 360.0 - diskrim.PointZeta[diskrim.points];
      diskrim.PointH[1] = ( zetaN * diskrim.PointH[2]
	                  + zeta2 * diskrim.PointH[diskrim.points])
                                   / (zeta2 + zetaN);
      diskrim.PointT[1] = ( zetaN * diskrim.PointT[2]
                          + zeta2 * diskrim.PointT[diskrim.points])
                                   / (zeta2 + zetaN);
   }

   diskrim.Hmax = diskrim.PointH[1];
   diskrim.Hmin = diskrim.PointH[1];
   diskrim.ZetaHmax = diskrim.PointZeta[1];
   diskrim.Tmin = diskrim.PointT[1];
   diskrim.Tmax = diskrim.PointT[1];
   diskrim.ZetaTmax = diskrim.PointZeta[1];
   for( i = 1; i <= diskrim.points; i++) {
      if( diskrim.PointH[i] > diskrim.Hmax ) {
	 diskrim.Hmax = diskrim.PointH[i];
         diskrim.ZetaHmax = diskrim.PointZeta[i];
      }
      if( diskrim.PointH[i] < diskrim.Hmin ) {
	 diskrim.Hmin = diskrim.PointH[i];
      }
      if( diskrim.PointT[i] > diskrim.Tmax ) {
	 diskrim.Tmax = diskrim.PointT[i];
         diskrim.ZetaTmax = diskrim.PointZeta[i];
      }
      if( diskrim.PointT[i] < diskrim.Tmin ) {
	 diskrim.Tmin = diskrim.PointT[i];
      }
   }

   return;
}


void MakeDiskTorus( void )
/******************************************************
*
*   This function is invoked if the disk torus is ON and the torus
*   is specified by a set of discrete points around the disk.  
*   It does the following:
*      1)  Sort the torus points in order of increasing zeta.
*      2)  Check for duplicate points.
*      3)  Create a point a zeta = 0 if one does not already
*            exist by linearly interpolating in zeta.
*      4)  Finds maximum and minimum values for H and T.
*
*   Note:  The function assumes that zeta is given in degrees.
*
*************************************************************/
{
   long i, j;
   double dummy, zeta2, zetaN;

   if( disktorus.points <= 0 )
      Quit("SortTorusPoints:  disktorus.points le 0.");
   else if( disktorus.points > 1 ) {
      for( i = 1; i < disktorus.points; i++ ) {
         for( j = i+1; j <= disktorus.points; j++) {
            if( disktorus.PointZeta[i] > disktorus.PointZeta[j] ) {
 	       dummy = disktorus.PointZeta[i];
               disktorus.PointZeta[i] = disktorus.PointZeta[j];
               disktorus.PointZeta[j] = dummy;
               dummy = disktorus.PointH[i];
               disktorus.PointH[i] = disktorus.PointH[j];
               disktorus.PointH[j] = dummy;
               dummy = disktorus.PointT[i];
               disktorus.PointT[i] = disktorus.PointT[j];
               disktorus.PointT[j] = dummy;
            }
         }
      }
   }
   if( disktorus.points > 1 ) {
      for( i = 1; i <= disktorus.points; i++) {
         if( disktorus.PointZeta[i] == disktorus.PointZeta[i+1] )
            Quit("DISKTORUSPARS: At least two points have the same PointZeta.");
      }
   }

   if( disktorus.PointZeta[1] != 0.0 ) {
      disktorus.points += 1;
      for( i = disktorus.points; i >= 2; i--) {
         disktorus.PointZeta[i] = disktorus.PointZeta[i-1];
         disktorus.PointH[i]    = disktorus.PointH[i-1];
         disktorus.PointT[i]    = disktorus.PointT[i-1];
      }
      disktorus.PointZeta[1] = 0.0;
      zeta2 = disktorus.PointZeta[2];
      zetaN = 360.0 - disktorus.PointZeta[disktorus.points];
      disktorus.PointH[1] = ( zetaN * disktorus.PointH[2]
	                  + zeta2 * disktorus.PointH[disktorus.points])
                                   / (zeta2 + zetaN);
      disktorus.PointT[1] = ( zetaN * disktorus.PointT[2]
                          + zeta2 * disktorus.PointT[disktorus.points])
                                   / (zeta2 + zetaN);
   }

   disktorus.Hmax = disktorus.PointH[1];
   disktorus.Hmin = disktorus.PointH[1];
   disktorus.ZetaHmax = disktorus.PointZeta[1];
   disktorus.Tmin = disktorus.PointT[1];
   disktorus.Tmax = disktorus.PointT[1];
   disktorus.ZetaTmax = disktorus.PointZeta[1];
   for( i = 1; i <= disktorus.points; i++) {
      if( disktorus.PointH[i] > disktorus.Hmax ) {
	 disktorus.Hmax = disktorus.PointH[i];
         disktorus.ZetaHmax = disktorus.PointZeta[i];
      }
      if( disktorus.PointH[i] < disktorus.Hmin ) {
	 disktorus.Hmin = disktorus.PointH[i];
      }
      if( disktorus.PointT[i] > disktorus.Tmax ) {
	 disktorus.Tmax = disktorus.PointT[i];
         disktorus.ZetaTmax = disktorus.PointZeta[i];
      }
      if( disktorus.PointT[i] < disktorus.Tmin ) {
	 disktorus.Tmin = disktorus.PointT[i];
      }
   }

   return;
}

