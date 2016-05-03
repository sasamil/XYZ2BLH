/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   sasa.milenkovic.xyz@gmail.com                                         *
 *                                                                         *
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

#ifndef _lmjjsAux
#define _lmjjsAux

typedef double _double_;
//-------------------------------------------------------------------------------------------

  const _double_ tolerance  =  1.e-12;
  const _double_ mtolerance = -tolerance;
  const _double_ dtolerance =  tolerance * tolerance;

  const _double_ PI = 3.141592653589793238463;
  const _double_ TWOPI = PI * 2.0;
  const _double_ HALFPI = PI * 0.5;
  const _double_ ro = 3600.0 * 180.0 / PI; // rad2sec
  const _double_ ONETHIRD = 1./3.;

  // parameters I often need. constants
  const _double_  a_bessel	   = 6377397.155;
  const _double_  a_bessel_2	 = a_bessel*a_bessel;
  const _double_  b_bessel	   = 6356078.962818189;
  const _double_  b_bessel_2	 = b_bessel*b_bessel;
  const _double_  f_bessel	   = 1.0L - b_bessel/a_bessel; //1.0L/299.1528128L;
  const _double_  f_bessel2	   = a_bessel/b_bessel - 1.0;
  const _double_  j_e2_bessel  = b_bessel_2/a_bessel_2;
  const _double_  e2_bessel	   = 1.0L - j_e2_bessel; //0.0816968312225269L;
  const _double_  e2_bessel2	 = (a_bessel_2 - b_bessel_2) / b_bessel;

  const _double_  a_wgs	    = 6378137.0;
  const _double_  a_wgs_2   = a_wgs*a_wgs;
  const _double_  b_wgs	    = 6356752.314245179;
  const _double_  b_wgs_2   = b_wgs*b_wgs;
  const _double_  f_wgs	    = (a_wgs - b_wgs) / a_wgs; //1.0L/298.257223563L;
  const _double_  f_wgs2	  = (a_wgs - b_wgs) / b_wgs;
  const _double_  j_e2_wgs  = b_wgs_2/a_wgs_2;
  const _double_  e2_wgs	  = 1.0 - j_e2_wgs; //0.08181919084262149L;
  const _double_  e2_wgs2   = (a_wgs_2 - b_wgs_2) / b_wgs;

  struct XYZtriplet {_double_ x; _double_ y; _double_ z;};

//-------------------------------------------------------------------------------------------
  unsigned int solveP3(_double_*, _double_, _double_, _double_);
  unsigned int solveP4(_double_*, _double_*, _double_, _double_, _double_, _double_);

  XYZtriplet BLh2XYZ(const _double_, const _double_, const XYZtriplet&);
  XYZtriplet geocentric_to_geodetic(const _double_, const _double_, const XYZtriplet&);
  XYZtriplet geocentric_to_geodetic_bowring (const _double_, const _double_, const XYZtriplet&);
  XYZtriplet bowring2(const _double_, const _double_, const XYZtriplet&);
  XYZtriplet bowring3(const _double_, const _double_, const XYZtriplet&);
  XYZtriplet nr2(const _double_ a, const _double_ b, const XYZtriplet& p);
  XYZtriplet sq2(const _double_ a, const _double_ b, const XYZtriplet& p);
  XYZtriplet direct_solution(const _double_ a, const _double_ b, const XYZtriplet& p);
  XYZtriplet direct_solution_opt(const _double_ a, const _double_ b, const XYZtriplet& p);

  //---------------------------------------------------------------------------

#endif
