#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;
using namespace LHAPDF;
const double conv = 1E-27;	//conversion from GeV^-2 to cm^2

/* The PDF Set and Set Member being used */
string setname = "CT10";
int setmember = 1;
/* Call PDF */
PDF* pdf = mkPDF(setname, setmember);

/* structure to hold all parameters that may change */
struct arg1_params{
	double Enu;	//neutrino energy
};


/* function that returns the differential cross-section*/
double func(double *x, size_t dim, void * p){
	
	/*parameters*/
	struct arg1_params * fp = (struct arg1_params*) p;
	double Enu = fp->Enu;
	double pi = M_PI;
	const double GF = 1.16632E-5;	//GeV^-2
	const double M = 0.9;		//isoscalar effective mass in GeV
	const double Mw = 80.0;		//GeV
	double Q2 = 2.0 * M * Enu * x[0] * x[1];

	/* Make a PDF for each quark/anti-quark type. PIDs are:
	down	up	strange	charm	beauty	truth
	1	2	3	4	5	6 */	
	double u = pdf->xfxQ2(2, x[0], Q2);
	double d = pdf->xfxQ2(1, x[0], Q2);
	double s = pdf->xfxQ2(3, x[0], Q2);
	double cbar = pdf->xfxQ2(-4, x[0], Q2);
	double b = pdf->xfxQ2(5, x[0], Q2);
	double ubar = pdf->xfxQ2(-2, x[0], Q2);
	double dbar = pdf->xfxQ2(-1, x[0], Q2);

	/*distributions*/
	double xq = (u-ubar + d-dbar)/2.0 + (ubar + dbar)/2.0 + s + b;
	double xqbar = (ubar + dbar)/2.0 + cbar;

	/*differential cross section*/
	double eta = ( 2.0*M*Enu/pi ) * pow(GF, 2) * pow( ( Mw*Mw/(Q2 + Mw*Mw) ), 2 );
	double f = eta* ( xq + xqbar * pow( (1.0 - x[1]), 2) ) ;
	return f;
}




int main(){

	/*use these to time the program*/
	clock_t ti, tf;
	ti = clock();
	
	cout << "initializing... " << '\n';

	/* relevant variables */
	double Sigmacc = 0.0;
	double Err = 0.0;
	size_t calls = 10000000;	//# function calls for VEGAS
	
	/* initial probe energy in GeV*/
	double Eneutrino = 10.0;

	/* integration bounds for {x, y} */
	const double xlow = pdf->xMin();
	const double ylow = 1E-8;
	const double xup = pdf->xMax();
	const double yup = 1.0;
	double lower[2] = {xlow, ylow};	
	double upper[2] = {xup, yup};

	/*initialize PDF stuff */
	setVerbosity(SILENT);
	initPDFSet(setname, setmember);	

	/*set up VEGAS stuff and random number generator*/
	gsl_monte_function F;
	gsl_monte_function TF;
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	gsl_rng *r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	/*initialize functions for integration */
	F.f = &func;
	F.dim = 2;

	/* open files to write */
	ofstream myfile1;
	ofstream myfile2;
	myfile1.open("NuN1.dat");
	myfile2.open("NuN2.dat");
	int i = 0; //dummy counter

	/* START INTEGRATING, keep cross-section as function of probe energy */
	while (Eneutrino <= 1E15){
		
		/*integration workspace*/
		gsl_monte_vegas_state * s =  gsl_monte_vegas_alloc(2);
		struct arg1_params P = {Eneutrino};
		F.params = &P;
		
		/* check status */
		if (i == 0){
			cout << "working..." << '\n';
			i++;
		}
		
		/*integration over Bjorken variables*/
		gsl_monte_vegas_integrate(&F, lower, upper, 2, calls, r, s, &Sigmacc, &Err);
		/* write to files */
		myfile1 <<(Eneutrino) << " " << log10(2E11 * Sigmacc/Eneutrino) << '\n';
		myfile2 << (Eneutrino) << " " << log10(conv * Sigmacc) << '\n';

		/* set probe energy, perpare to start again */
		Eneutrino = Eneutrino * 10.0;
		gsl_monte_vegas_free(s);
	}

	/* close files, clear memory, finish timing */
	myfile1.close();
	myfile2.close();
	gsl_rng_free(r);
	tf = clock();
	double tdiff = (double)tf - (double)ti;
	double sec = tdiff / CLOCKS_PER_SEC;
	double mins = sec/60.0;
	if (sec < 60){
		cout << '\n' << "Time: " << sec << "seconds" << '\n' << '\n';
	}

	else {
		cout << '\n' << "Time: " << mins << "minutes" << '\n' << '\n';
	}

	return 0;

}

