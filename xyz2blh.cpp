/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   sasa.milenkovic.xyz@gmail.com                                         *
 *   									   *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *   ( http://www.gnu.org/licenses/gpl-3.0.en.html )                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cmath>
#include <iostream>

#include "xyz2blh.h" //all includes are here


//---------------------------------------------------------------------------
//input: fi, lambda (in radians) and h (m)
//output: geocentric coordinates (m)
XYZtriplet BLh2XYZ(const _double_ a, const _double_ b, const XYZtriplet& blh)
{
	_double_ fi = blh.x;
	_double_ lambda = blh.y;
	_double_ h = blh.z;

	_double_ e_2 = 1.0 - b*b/a/a;
	_double_ sn = sin(fi);
	_double_ cs = cos(fi);
	_double_ N = a/sqrt(1.0 - e_2*sn*sn);

	_double_ X = (N + h)		* cs*cos(lambda);
	_double_ Y = (N + h)		* cs*sin(lambda);
	_double_ Z = (N*(1-e_2) + h)	* sn;

	XYZtriplet retval;
	retval.x = X;
	retval.y = Y;
	retval.z = Z;

	return retval;
}

//--------------------------------------------------------------------------------
//input: geocentric coordinates (m)
//output: fi, lambda (in radians) and h (m)
XYZtriplet geocentric_to_geodetic(const _double_ a, const _double_ b, const XYZtriplet& p)
{
/* local defintions and variables */
/* end-criterium of loop, accuracy of sin(Latitude) */
   int maxiter = 30;
   _double_ genau = tolerance; //1.E-12;
   _double_ genau2 = (genau*genau);
   _double_ e_2 = 1.0 - b*b/a/a;

    bool At_Pole;    /* indicates location is in polar region */
    int iter;        /* # of continous iteration, max. 30 is always enough (s.a.) */
    _double_ P;        /* distance between semi-minor axis and location */
    _double_ RR;       /* distance between center and location */
    _double_ CT;       /* sin of geocentric latitude */
    _double_ ST;       /* cos of geocentric latitude */
    _double_ RX;
    _double_ RK;
    _double_ RN;       /* Earth radius at location */
    _double_ CPHI0;    /* cos of start or old geodetic latitude in iterations */
    _double_ SPHI0;    /* sin of start or old geodetic latitude in iterations */
    _double_ CPHI;     /* cos of searched geodetic latitude */
    _double_ SPHI;     /* sin of searched geodetic latitude */
    _double_ SDPHI;    /* end-criterium: addition-theorem of sin(Latitude(iter)-Latitude(iter-1)) */
	
    _double_ X = p.x;
    _double_ Y = p.y;
    _double_ Z = p.z ? p.z : 0.0;   //Z value not always supplied
    _double_ Longitude;
    _double_ Latitude;
    _double_ Height;
    XYZtriplet retval;

    At_Pole = false;
    P = sqrt(X*X+Y*Y);
    RR = sqrt(X*X+Y*Y+Z*Z);

/*  special cases for latitude and longitude */
    if (P/a < genau) 
    {

/*  special case, if P=0. (X=0., Y=0.) */
        At_Pole = true;
        Longitude = 0.0;

/*  if (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis
 *  of ellipsoid (=center of mass), Latitude becomes PI/2 */
        if (RR/a < genau) {
        	Latitude = HALFPI;
        	Height   = -b;

		retval.x = Latitude;
		retval.y = Longitude;
		retval.z = Height;
        	return retval;
        }
    }
    else
    {
/*  ellipsoidal (geodetic) longitude
 *  interval: -PI < Longitude <= +PI */
        Longitude=atan2(Y,X);
    }

/* --------------------------------------------------------------
 * Following iterative algorithm was developped by
 * "Institut für Erdmessung", University of Hannover, July 1988.
 * Internet: www.ife.uni-hannover.de
 * Iterative computation of CPHI,SPHI and Height.
 * Iteration of CPHI and SPHI to 10**-12 radian resp.
 * 2*10**-7 arcsec.
 * --------------------------------------------------------------
 */
    CT = Z/RR;
    ST = P/RR;
    RX = 1.0/sqrt(1.0-e_2*(2.0-e_2)*ST*ST);
    CPHI0 = ST*(1.0-e_2)*RX;
    SPHI0 = CT*RX;
    iter = 0;

/* loop to find sin(Latitude) resp. Latitude
 * until |sin(Latitude(iter)-Latitude(iter-1))| < genau */
    do
    {
        iter++;
        RN = a/sqrt(1.0-e_2*SPHI0*SPHI0);

/*  ellipsoidal (geodetic) height */
        Height = P*CPHI0+Z*SPHI0-RN*(1.0-e_2*SPHI0*SPHI0);

        RK = e_2*RN/(RN+Height);
        RX = 1.0/sqrt(1.0-RK*(2.0-RK)*ST*ST);
        CPHI = ST*(1.0-RK)*RX;
        SPHI = CT*RX;
        SDPHI = SPHI*CPHI0-CPHI*SPHI0;
        CPHI0 = CPHI;
        SPHI0 = SPHI;
    }
    while (SDPHI*SDPHI > genau2 && iter < maxiter);

/*      ellipsoidal (geodetic) latitude */
    Latitude=atan(SPHI/fabs(CPHI));

	retval.x = Latitude;
	retval.y = Longitude;
	retval.z = Height;
    return retval;
  } // geocentric_to_geodetic()

//---------------------------------------------------------------------------
// equation to solve:  x^3 + a*x^2 + b*x + c
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(_double_* x, _double_ a, _double_ b, _double_ c)
{
	_double_ a2 = a*a;
    	_double_ q  = (a2 - 3.*b)/9.;
	_double_ r  = (a*(2.*a2-9.*b) + 27.*c)/54.;
    	_double_ r2 = r*r;
	_double_ q3 = q*q*q;
	_double_ A,B;
	
    	if(r2<q3) 
    	{
        	_double_ t=r/sqrt(q3);
		if( t<-1.) t=-1.;
		if( t> 1.) t= 1.;
        	t=acos(t);
        	a/=3.; q=-2.*sqrt(q);
        	x[0]=q*cos(t/3.)-a;
        	x[1]=q*cos((t+TWOPI)/3.)-a;
        	x[2]=q*cos((t-TWOPI)/3.)-a;
        	return(3);
    	} 
    	
    	else 
    	{
        	//A =-pow(fabs(r)+sqrt(r2-q3),1./3.);
        	A = -exp(log(fabs(r)+sqrt(r2-q3)) * ONETHIRD);
		if( r<0. ) A=-A;
		B = (0.==A ? 0. : q/A);

		a/=3.;
		x[0] =(A+B)-a;
        	x[1] =-0.5*(A+B)-a;
        	x[2] = 0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<tolerance) { x[2]=x[1]; return(2); }
        	return(1);
	 }
}// solveP3

//---------------------------------------------------------------------------
// equation to solve:  x^4 + a*x^3 + b*x62 + c*x + d
// re - array of size 4, im - array of size 2
// returns number of
// In case 4 real roots:  re[0], re[1], re[2], re[3]                                              return 4
//         2 real roots:  re[0], re[1], re[2] = re[3]  or   x[0] ± i*x[1], x[2], x[3]             return 2
//         0 real roots : x[0] ± i*x[1], x[2] ± i*x[3]                                            return 0
// In case x[1] or x[3] are imaginnary -
unsigned int solveP4(_double_ *re, _double_ *im, _double_ a, _double_ b, _double_ c, _double_ d)
{
    	_double_ q1, q2, p1, p2, D, sqd, y;
    	_double_ a3 = -b;
    	_double_ b3 =  a*c -4.*d;
	_double_ c3 = -a*a*d - c*c + 4.*b*d;
	_double_ x3[3];

	unsigned int iRetval=4, iZeroes = solveP3(x3, a3, b3, c3);

	y = x3[0];
	if(iZeroes != 1.)
	{
		if(fabs(x3[1]) > fabs(y)) y = x3[1];
		if(fabs(x3[2]) > fabs(y)) y = x3[2];
	}

	D = y*y - 4.*d;
	if(D<tolerance /*fabs(D) < eps*/)
	{
		q1 = q2 = y * 0.5;

		D = a*a - 4.*(b-y);
		if(D<tolerance /*fabs(D) < eps*/)
			p1 = p2 = a * 0.5;
		else
		{
			sqd = sqrt(D);
			p1 = (a + sqd/*sqrt(D)*/) * 0.5;
			p2 = (a - sqd/*sqrt(D)*/) * 0.5;
		}
	}
	else
	{
		sqd = sqrt(D);
        	q1 = (y + sqd/*sqrt(D)*/) * 0.5;
		q2 = (y - sqd/*sqrt(D)*/) * 0.5;

		p1 = (a*q1-c)/(q1-q2);
		p2 = (c-a*q2)/(q1-q2);
	}

	im[0] = im[1] = 0.0;

	// x^2 + p1*x + q1 = 0
	D = p1*p1 - 4.*q1;
	if(D < 0.0)
	{
		iRetval -= 2.;
        	re[0] = re[1] = -p1 * 0.5;
        	im[0] = sqrt(-D) * 0.5;
	}
	else
	{
		sqd = sqrt(D);
        	re[0] = (-p1 + sqd/*sqrt(D)*/) * 0.5;
        	re[1] = (-p1 - sqd/*sqrt(D)*/) * 0.5;
	}

	// x^2 + p2*x + q2 = 0
	D = p2*p2 - 4.*q2;
	if(D < 0.0)
	{
		iRetval -= 2.;
        	re[2] = re[3] = -p2 * 0.5;
        	im[1] = sqrt(-D) * 0.5;
	}
	else
	{
        	sqd = sqrt(D);
        	re[2] = (-p2 + sqd/*sqrt(D)*/) * 0.5;
        	re[3] = (-p2 - sqd/*sqrt(D)*/) * 0.5;
	}

	return iRetval;
}// solveP4

//---------------------------------------------------------------------------
// The method used here is derived from 'An Improved Algorithm for
// Geocentric to Geodetic Coordinate Conversion, by Ralph Toms, Feb 1996
// Note: Variable names follow the notation used in Toms, Feb 1996
#define PI_OVER_2  (3.14159265358979323e0 / 2.0e0)
#define FALSE      0
#define TRUE       1
#define COS_67P5   0.38268343236508977  /* cosine of 67.5 degrees */
#define AD_C       1.0026000            /* Toms region 1 constant */

XYZtriplet geocentric_to_geodetic_bowring (const _double_ Geocent_a, const _double_ Geocent_b, const XYZtriplet& p)
{ /* BEGIN Convert_Geocentric_To_Geodetic */
    _double_ W;        /* distance from Z axis */
    _double_ W2;       /* square of distance from Z axis */
    _double_ T0;       /* initial estimate of vertical component */
    _double_ T1;       /* corrected estimate of vertical component */
    _double_ S0;       /* initial estimate of horizontal component */
    _double_ S1;       /* corrected estimate of horizontal component */
    _double_ Sin_B0;   /* sin(B0), B0 is estimate of Bowring aux variable */
    _double_ Sin3_B0;  /* cube of sin(B0) */
    _double_ Cos_B0;   /* cos(B0) */
    _double_ Sin_p1;   /* sin(phi1), phi1 is estimated latitude */
    _double_ Cos_p1;   /* cos(phi1) */
    _double_ Rn;       /* Earth radius at location */
    _double_ Sum;      /* numerator of cos(phi1) */
    int At_Pole;     /* indicates location is in polar region */

    /* sm */
    _double_ X = p.x;
    _double_ Y = p.y;
    _double_ Z = p.z;
    _double_ Geocent_a2 = Geocent_a * Geocent_a;
    _double_ Geocent_b2 = Geocent_b * Geocent_b;
    _double_ Geocent_e2 = (Geocent_a2 - Geocent_b2) / Geocent_a2;
    _double_ Geocent_ep2 = (Geocent_a2 - Geocent_b2) / Geocent_b2;

    XYZtriplet blh;

    At_Pole = FALSE;
    if (X != 0.0)
    {
        //*Longitude = atan2(Y,X);
        blh.y = atan2(Y,X);
    }
    else
    {
        if (Y > 0)
        {
            //*Longitude = PI_OVER_2;
            blh.y = PI_OVER_2;
        }
        else if (Y < 0)
        {
            //*Longitude = -PI_OVER_2;
            blh.y = -PI_OVER_2;
        }
        else
        {
            At_Pole = TRUE;
            //*Longitude = 0.0;
            blh.y = 0.0;

            if (Z > 0.0)
            {  /* north pole */
                //*Latitude = PI_OVER_2;
                blh.x = PI_OVER_2;
            }
            else if (Z < 0.0)
            {  /* south pole */
                //*Latitude = -PI_OVER_2;
                blh.x = -PI_OVER_2;
            }
            else
            {  /* center of earth */
                //*Latitude = PI_OVER_2;
                //*Height = -Geocent_b;
                blh.x = PI_OVER_2;
                blh.z = -Geocent_b;
                return blh;
            }
        }
    }
    
    W2 = X*X + Y*Y;
    W = sqrt(W2);
    T0 = Z * AD_C;
    S0 = sqrt(T0 * T0 + W2);
    Sin_B0 = T0 / S0;
    Cos_B0 = W / S0;
    Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;
    T1 = Z + Geocent_b * Geocent_ep2 * Sin3_B0;
    Sum = W - Geocent_a * Geocent_e2 * Cos_B0 * Cos_B0 * Cos_B0;
    S1 = sqrt(T1*T1 + Sum * Sum);
    Sin_p1 = T1 / S1;
    Cos_p1 = Sum / S1;
    Rn = Geocent_a / sqrt(1.0 - Geocent_e2 * Sin_p1 * Sin_p1);

    blh.x = atan(Sin_p1 / Cos_p1);

    if (Cos_p1 >= COS_67P5)
    {
        //*Height = W / Cos_p1 - Rn;
        blh.z = W / Cos_p1 - Rn;
    }
    else if (Cos_p1 <= -COS_67P5)
    {
        //*Height = W / -Cos_p1 - Rn;
        blh.z = W / -Cos_p1 - Rn;
    }
    else
    {
        //*Height = Z / Sin_p1 + Rn * (Geocent_e2 - 1.0);
        blh.z = Z / Sin_p1 + Rn * (Geocent_e2 - 1.0);
    }
    if (At_Pole != FALSE)
    {
        //*Latitude = atan(Sin_p1 / Cos_p1);
        blh.z = atan(Sin_p1 / Cos_p1);
    }

    return blh;
} // geocentric_to_geodetic_bowring

//---------------------------------------------------------------------------
// Bowring: /home/sasamil/Geodesy/BLH2XYZ/COMPARISON OF DIFFERENT ALGORITHMS.pdf  or
// /home/sasamil/Geodesy/Books/GPS_Hofmann-Wellenhof (p.232)   or
// /home/sasamil/Geodesy/BLH2XYZ/Comparision2!! (p.8)
XYZtriplet bowring2(const _double_ a, const _double_ b, const XYZtriplet& p)
{
	_double_ f2  = a/b;
	_double_ epr2 = f2*f2 - 1.0;
	_double_ e2 = (a*a - b*b)/(a*a);
	_double_ X = p.x;
	_double_ Y = p.y;
	_double_ Z = p.z;
	_double_ P = sqrt(X*X + Y*Y);

	_double_ tgtheta, cstheta, sntheta, snfi, N, tgfi, tgfiold, delta;
	XYZtriplet blh;

    if(0.0==X && 0.0==Y)
    {
        blh.x = HALFPI;
        blh.y = 0.0;
        blh.z = Z - b;
        return blh;
    }
    else if(0.0==Z)
    {
        blh.x = 0.0;
        blh.y = atan(Y/X);
        blh.z = P - a;
        return blh;
    }

    else
    {
        tgfi = f2*f2*Z/P;
        do
        {
                tgfiold = tgfi;
                tgtheta = tgfi/f2;
                //_double_ theta = atan(f2*Z/P);
                cstheta = 1./sqrt(1. + tgtheta*tgtheta); //cos(theta);
                sntheta = tgtheta*cstheta; //tg/sqrt(1. + tg*tg); //sin(theta);

                tgfi = (Z + epr2*b*sntheta*sntheta*sntheta)/(P - e2*a*cstheta*cstheta*cstheta);
                //fi = atan(tg2);
                //snfi = sin(fi);
                delta = tgfi-tgfiold;
        }
        while(delta*delta > dtolerance);

        snfi = tgfi / sqrt(1. + tgfi*tgfi);
        N = a / sqrt(1. - e2*snfi*snfi);

        blh.x = atan(tgfi);  //E = atan( (1-f_bessel)*tan(B) );
        blh.y = atan(Y/X);
        blh.z = P*tgfi/snfi - N;

        return blh;
    }
} // bowring


//---------------------------------------------------------------------------
// Bowring: http://gis.stackexchange.com/questions/28446/..
// ..computational-most-efficient-way-to-convert-cartesian-to-geodetic-coordinates
XYZtriplet bowring3(const _double_ a, const _double_ b, const XYZtriplet& pin)
{
    _double_ X = pin.x;
    _double_ Y = pin.y;
    _double_ Z = pin.z;
    _double_ a2 = a*a;
    _double_ f2  = a/b;
    _double_ eb2 = f2*f2 - 1.0;
    _double_ e2 = (a2 - b*b)/a2;
    _double_ d  = eb2*b;
    //_double_ ome2 = 1.0 - e2;
    _double_ p, tu, tu2, su3, cu, cu3, tp, cp, sp, delta, tpold;

    XYZtriplet blh;

    if(0.0==X && 0.0==Y)
    {
        blh.x = HALFPI;
        blh.y = 0.0;
        blh.z = Z - b;
        return blh;
    }
    else if(0.0==Z)
    {
        blh.x = 0.0;
        blh.y = atan(Y/X);
        p = sqrt(X*X + Y*Y);
        blh.z = p - a;
        return blh;
    }
    else
    {
        p = sqrt(X*X + Y*Y);
        tp = f2*f2*Z/p;
        do
        {
            tpold = tp;
            tu  = tp/f2; //b*Z*(1.0 + d/r)/(a*p);
            tu2 = tu*tu;
            cu  = (1.0/sqrt(1.0 + tu2));
            cu3 = cu*cu*cu;
            su3 = cu3*tu2*tu;
            tp  = (Z + d*su3)/(p - e2*a*cu3);
            delta = tp-tpold;
        }
        while(delta*delta > dtolerance);

        cp  = 1.0/sqrt(1.0 + tp*tp);
        sp  = cp*tp;

        blh.x = atan(tp);  //E = atan( (1-f_bessel)*tan(B) );
        blh.y = atan(Y/X);
        blh.z = p*cp + Z*sp - a*sqrt(1.0 - e2*sp*sp);

        return blh;
    }
} // bowring2


//---------------------------------------------------------------------------
// the basis is: (1+f')*X*tgE - b*e'^2*sinE - Y = 0   
// Solving system y = x^4 + A*t^3 + B*t - 1 by Newton-Raphson ( t=tg(E/2) )
XYZtriplet nr2(const _double_ a, const _double_ b, const XYZtriplet& p)
{
    //_double_ tolerance = 1.E-12;
    _double_ step, x1, x1p2, y1, y1pr,
        f2, e22, f2R, X, Y, Z, R, A, B, C, tgEpocetno,
        tg, cs;

    XYZtriplet blh;

    X = p.x;
    Y = p.y;
    Z = p.z;
    R = sqrt(X*X + Y*Y);

    if(0.0==X && 0.0==Y)
    {
        blh.x = HALFPI;
        blh.y = 0.0;
        blh.z = Z - b;
        return blh;
    }
    else if(0.0==Z)
    {
        blh.x = 0.0;
        blh.y = atan(Y/X);
        blh.z = R - a;
        return blh;
    }

    else
    {
        f2  = a/b;
        e22 = a*f2 - b;
        f2R = f2*R;

        A = 2.*(f2R + e22) / Z;
        B = 2.*(f2R - e22) / Z;
        C = 3.*A;

        tgEpocetno = f2*Z/R;
        x1 = tgEpocetno / (1.+sqrt(1.+tgEpocetno*tgEpocetno));  // x=t=tan(E/2);
        x1p2 = x1*x1; // Die erste Iteration ist ausgeschlossen, weil Funktion könnte negativ sein (ausschließlich).

        // Visual Studio
        //y1 = x1*(B + x1p2*(A + x1)) - 1.;
        //y1pr = x1p2*(C + 4.*x1) + B;
        //step = y1/y1pr;

        // gcc - linux
        y1 = x1;
        y1pr = x1;
        step = --((((y1+=A)*=x1p2)+=B)*=x1) / ((((y1pr*=4.)+=C)*=x1p2)+=B);
        
        x1 -= step;

        do
        {
            x1p2 = x1*x1; // Die erste Iteration ist ausgeschlossen, weil Funktion könnte negativ sein (ausschließlich).

            // Visual Studio
            //y1 = x1*(B + x1p2*(A + x1)) - 1.;
            //y1pr = x1p2*(C + 4.*x1) + B;
            //step = y1/y1pr;

            // gcc - linux
            y1 = y1pr = x1;
            step = --((((y1+=A)*=x1p2)+=B)*=x1) / ((((y1pr*=4.)+=C)*=x1p2)+=B);
            
            x1 -= step;
        }
        while (step > tolerance);

        tg = (2.*x1)/(1.-x1*x1) * f2;
        cs = 1. / sqrt(1.+tg*tg);
        //t = (1.-t)/(1.+t); // t = tg(pi/4 - E/2) -->  t = tg(E/2)
        blh.x = atan(tg);
        blh.y = atan(Y/X);
        blh.z = (R - a*(1.-x1)/(1.+x1))*cs + (Z-b)*tg*cs;

        return blh;
    }
} // nr2

//---------------------------------------------------------------------------
// Solving system y = x^4 + A*t^3 + B*t - 1 by quadratic interpolation
// t = tg(E/2)
XYZtriplet sq2(const _double_ a, const _double_ b, const XYZtriplet& p)
{
	_double_ f2, e22, X, Y, Z, R, A, B,
            x1, x1p2, y1, x2, x2p2, y2, oldx2, oldy,
            tg, cs;


    	XYZtriplet blh;

	f2  = a/b;
	e22 = a*f2 - b;
	X = p.x;
	Y = p.y;
	Z = p.z;
	R = sqrt(X*X + Y*Y);

	A = 2.*(f2*R + e22) / Z;
	B = 2.*(f2*R - e22) / Z;

	x1 = Z/R * f2;  x1 = x1 / (1. + sqrt(1.+x1*x1));
	x1p2 = x1*x1;   y1 = x1*(B + x1p2*(A + x1)) - 1.;

	x2 = x1 - y1/(x1p2*(3.*A + 4.*x1) + B);
    	x2p2 = x2*x2;   y2 = x2*(B + x2p2*(A + x2)) - 1.;

	oldx2=x1p2, oldy=y1;
	do
	{
		x1p2 = (x2p2*y1 - x1p2*y2) / (y1-y2);
		x1   = sqrt(x1p2);
		y1=x1*(B + x1p2*(A + x1)) - 1.;

		x2p2 = oldx2, y2 = oldy;
		oldx2 = x1p2,  oldy = y1;
    	}
    	while( (y1>tolerance) | (y1<mtolerance) );

    	tg = (2.*x1)/(1.-x1*x1) * f2;
    	cs = 1. / sqrt(1.+tg*tg);
    	//t = (1.-t)/(1.+t); // t = tg(pi/4 - E/2) -->  t = tg(E/2)
    	blh.x = atan(tg);
    	blh.y = atan(Y/X);
    	blh.z = (R - a*(1.-x1)/(1.+x1))*cs + (Z-b)*tg*cs;

	return blh;
} // sq2

//---------------------------------------------------------------------------
// the basis is: (1+f')*X*tgE - b*e'^2*sinE - Y = 0   
// Solving system y = x^4 + A*t^3 + B*t^2 + C*t + D by direct method  ( t=tg(E/2) )
XYZtriplet direct_solution(const _double_ a, const _double_ b, const XYZtriplet& p)
{
	//_double_ f2  = a/b;
	//_double_ e22 = a*f2 - b;
	_double_ X = p.x;
	_double_ Y = p.y;
	_double_ Z = p.z;
	_double_ R = sqrt(X*X + Y*Y);

	//_double_ A = 2*(f2*R + e22) / Z;
	//_double_ B = 2*(f2*R - e22) / Z;

    	_double_ f  = b/a;
    	_double_ sqd = -f * Z/R;
    	_double_ pom = (f*b-a)/R;

    	_double_ A = 2. * sqd;
    	_double_ C = A;
    	_double_ D = sqd*sqd;
    	_double_ B = 1 - pom*pom + D;

	_double_ re[4], im[2];

	solveP4(re, im, A, B, C, D);

	_double_ t;
	XYZtriplet blh;
	if(0.0==im[0])
        	t = re[0]>0.0 ? re[0] : re[1];
    	else
        	t = re[2]>0.0 ? re[2] : re[3];

    	///*
    	_double_ sn = t/sqrt(1. + t*t);
    	//_double_ cs = sn/x1;		x1=tg(E)

    	_double_ xe = a*sn/t;
    	_double_ ye = b*sn;

    	_double_ tgphi = t/f;
    	_double_ cosphi = 1.  / sqrt(1.+tgphi*tgphi);

    	blh.x = atan(tgphi);
    	blh.y = atan(Y/X);
    	blh.z = (R-xe)*cosphi + (Z-ye)*tgphi*cosphi;
    	//*/

	return blh;
} //direct_solution

//*---------------------------------------------------------------------------
/// the basis is: (1+f')*X*tgE - b*e'^2*sinE - Y = 0
// Solving system y = t^4 + A*t^3 + B*t - 1 by direct method
XYZtriplet direct_solution_opt(const double a, const double b, const XYZtriplet& p)
{
	double  f2, e22, X, Y, Z, R,
            aa, c, q1, q2, p1, p2,
            D, sqd, y, q, r, A, B,
            t, t2, tg, sn, xe, ye, tgphi, cosphi;

	f2  = a/b;
	e22 = a*f2 - b;
	X = p.x;
	Y = p.y;
	Z = p.z;
	R = sqrt(X*X + Y*Y);

	aa = 2*(f2*R + e22) / Z;
	c  = 2*(f2*R - e22) / Z;

    	q  = -(aa*c + 4.)/3.;
	r  =  (aa*aa - c*c)*0.5;
    	A  = -pow(fabs(r)+sqrt(r*r-q*q*q),1./3.);
    	//A = -exp(log(fabs(r)+sqrt(r*r-q*q*q)) * ONETHIRD);
    	if( r<0 ) A=-A;
    	B = (0==A ? 0 : q/A);
	y = A+B;

	D = y*y + 4.;
	sqd = sqrt(D);
    	q1 = (y + sqd) * 0.5;
    	q2 = (y - sqd) * 0.5;

    	p1 = (aa*q1-c)/(q1-q2);

	D = p1*p1 - 4*q1;
	if(D >= 0.0)
        	t = (-p1 + sqrt(D)) * 0.5;
	else
	{
        	p2 = (c-aa*q2)/(q1-q2);
        	t = (-p2 + sqrt(p2*p2 - 4*q2)) * 0.5; // D2 = p2*p2 - 4*q2;
    	}

    	XYZtriplet blh;
    	
	t2 = t*t;
	t *= 2;
	tg = t/(1-t2);
	sn = t/(1+t2);

	xe = a*sn/tg;
	ye = b*sn;

	tgphi = tg*f2;
	blh.x = atan(tgphi);			//E = atan( (1-f_bessel)*tan(B) );
	blh.y = atan(Y/X);
	cosphi = 1 / sqrt(1+tgphi*tgphi);
	blh.z = (R - xe)*cosphi + (Z-ye)*tgphi*cosphi;

	return blh;
}



