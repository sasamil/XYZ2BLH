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


#include <iostream>
#include <iomanip>	//setpreciosion

#include "xyz2blh.h"

//---------------------------------------------------------------------------
// This is how I test these things...
int main(int argc, char *argv[])
{
	XYZtriplet blh, xyz, blhret;
	std::cout << std::fixed << std::setprecision(12);

	///*
	/// Ртањ
	_double_ B = 43.7761 * PI/180.0;
	_double_ L = 21.8933 * PI/180.0;
	//_double_ H = 1570.0;
	_double_ H = 400000.0; // Dove
	//_double_ H = 36000000.0; // BeiDou-1
	//_double_ H = -5157000.0;

	blh.x = B;
	blh.y = L;
	blh.z = H;

	xyz = BLh2XYZ(a_bessel, b_bessel, blh);

	std::cout << "Bowring-proj test:" << std::endl;
	blhret = geocentric_to_geodetic_bowring(a_bessel, b_bessel, xyz);
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	std::cout << "Moritz-Heiskenen test:" << std::endl;
	blhret = geocentric_to_geodetic(a_bessel, b_bessel, xyz);
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	std::cout << "Bowring2 metoda test:" << std::endl;
	blhret = bowring2(a_bessel, b_bessel, xyz);
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	std::cout << "Bowring3 metoda test:" << std::endl;
	blhret = bowring3(a_bessel, b_bessel, xyz);
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	std::cout << "Newton-Raphson - nr2 test:" << std::endl;
	blhret = nr2(a_bessel, b_bessel, xyz);
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	std::cout << "Square interpolation -sq test" << std::endl;
	blhret = sq2(a_bessel, b_bessel, xyz);
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	blhret = direct_solution(a_bessel, b_bessel, xyz);
	std::cout << "Direct method test:" << std::endl;
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;

	blhret = direct_solution_opt(a_bessel, b_bessel, xyz);
	std::cout << "Optimized direct method test:" << std::endl;
	std::cout << blhret.x - B << std::endl;
	std::cout << blhret.y - L << std::endl;
	std::cout << blhret.z - H << std::endl;
	std::cout << std::endl;
	//*/

	//measuring time
	tm *local;
	time_t now, oldtime, elapsed;
	clock_t ticks;


	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of Moritz-Heiskenen:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = geocentric_to_geodetic(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of Moritz-Heiskenen:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of Bowring-proj:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = geocentric_to_geodetic_bowring(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of Bowring-proj:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of Bowring 2:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = bowring2(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of Bowring 2:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of Bowring 3:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;//+400000.0;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = bowring3(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of Bowring - 3:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of Newton-Raphson - nr2:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;//+400000.0;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = nr2(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of Newton-Raphson - nr2:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of quadratic interpolation - sq2:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;//+400000.0;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = sq2(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of quadratic interpolation - sq2:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of direct solution:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;//+400000.0;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = direct_solution(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of direct solution:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//--------------------------
	now = time(NULL);
	local = localtime(&now);
	std::cout << std::endl << "Before the millions executions of optimized direct solution:" << std::endl;
	std::cout << "Local Time ==> " << asctime(local) << std::endl;
	oldtime = now;
	ticks = clock();

	for(_double_ fi=0.1; fi<90.0; fi += 0.1)
		//for(_double_ h=0.0; h<10000.0; h+=1.0)
		for(_double_ h=395000; h<405000; h+=1.0)
		//for(_double_ h=35779000; h<35789000; h+=1.0)
		{
			blh.x = fi * PI/180.0;
			blh.z = h;//+400000.0;
			xyz = BLh2XYZ(a_bessel, b_bessel, blh);

			blhret = direct_solution_opt(a_bessel, b_bessel, xyz);
		}

	ticks = clock() - ticks;
	now = time(NULL);
	local = localtime(&now);
	elapsed = now-oldtime;
	std::cout << std::endl << "    " << "Ticks: " << ticks << std::endl;
	std::cout << "    " << "After the millions executions of optimized direct solution:" << std::endl;
	std::cout << "    " << "Local Time ==> " << asctime(local);// << endl;
	std::cout << "    " << "Elapsed Time: " << elapsed << " sec" << std::endl << std::endl;

	//std::cout << std::endl << "Das Program hat ein Scluß gemacht. Er ist erfolgreiche zu einem Ende gebracht!" << std::endl << std::endl;

	return 0;
}

