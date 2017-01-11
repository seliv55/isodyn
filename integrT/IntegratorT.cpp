/***************************************************************************
                               IntegratorT.cpp
                             -------------------
    written by           : Blake Ashby
    last updated         : Nov 15, 2002
    email                : bmashby@stanford.edu
 ***************************************************************************/


#include <fstream>
#include <iostream>
#include "IntegratorT.h"

// constructors

IntegratorT::IntegratorT(const int nin, double yin[], double xin, const double xendin,
	double dxin, int itolerin, double *rtolerin, double *atolerin, const int ioutin,
	double hin, double hmaxin, int nmaxin, double uroundin, double safein,
	double faclin, double facrin) :
	n(nin), y(yin), x(xin), xend(xendin), dx(dxin), itoler(itolerin),
	rtoler(rtolerin), atoler(atolerin), rtolerNULL(false), atolerNULL(false),
	iout(ioutin), h(hin), hmax(hmaxin),
	nmax(nmaxin), uround(uroundin), safe(safein), facl(faclin), facr(facrin),
	nfcn(0), nstep(0), naccpt(0), nrejct(0), xold(xin), hold(hin), xd(xin)
{
	// n, the dimension of the system
	if (n == UINT_MAX) {
		std::cout << "System too big, max. n = " << UINT_MAX - 1 << std::endl;
		throw -1;
	}
	
	// rtoler, the relative tolerance of the integration
	if (!rtoler) {
		itoler = 0;
		rtoler = new double;
		*rtoler = 1.0e-7;
		rtolerNULL = true;
	}
	
	// atoler, the absolute tolerance of the integration
	if (!atoler) {
		itoler = 0;
		atoler = new double;
		*atoler = 1.0e-7;
		atolerNULL = true;
	}

	// -------- maximal step size
	if (hmax == 0.0) hmax = xend - x;

	// -------- nmax--maximal number of steps
	if (nmax == 0) nmax = 100000;
	if (nmax <= 0) {
		std::cout << " wrong input, nmax = " << nmax << std::endl;
		throw -1;
	}
	
	// -------- uround--smallest number satisfying 1.0 + uround > 1.0
	if (uround == 0.0) uround = 1.0e-16;
	if ((uround <= 1.0e-19) || (uround >= 1.0)) {
		std::cout << " coefficients have 20 digits, uround = " << uround << std::endl;
		throw -1;
	}
	
	// --------- safe--safety factor in step size prediction
	if (safe == 0.0) safe = 0.9;
	if ((safe <= 0.001) || (safe >= 1.0)) {
		std::cout << " curious input for safety factor, safe = " << safe << std::endl;
		throw -1;
	}

}  // Constructor

// Destructor
IntegratorT::~IntegratorT()
{
	if (rtolerNULL) delete rtoler;
	if (atolerNULL) delete atoler;
}

// Function that controls the output of the results.
// Modify this routine according to your needs
int IntegratorT::SolutionOutput()
{
std::ofstream fout("kin");
	fout<<"t h6p= glgn= fbp= nt3= h6a= pep= pyr= oaa= mal= cit= coa= p5= e4= s7="<<std::endl;
//	std::cout << std::setiosflags(std::ios::showpoint);// | ios::fixed);

	if (naccpt == 0) xd = xold;

	while (xd <= x) 
		if (xold <= xd)  {fout <<xd<<" ";
/*			std::cout << "Step " << naccpt << ": t = " << std::setw(5) <<
				std::setprecision(2) << x << "  y = ";
			for (unsigned i = 0; i < n; i++)
				std::cout << std::setw(10) << std::setprecision(8)*/
				
for (unsigned i = 0; i < n; i++) fout << ContinuousOutput(i) << "  ";
			fout << std::endl;
			xd += dx;
		}
	

	return 0;

}  // SolutionOutput
