#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;
using namespace LHAPDF;
const double conv = 5E10;	//conversion from GeV^-2 to cm^2
const double conv2 = 7.77E-38/1.25736E-10;
/* The PDF Set and Set Member being used */
string setname = "CT10";
int setmember = 1;
/* Call PDF */
PDF* pdf = mkPDF(setname, setmember);

/* structure to hold all parameters that may change */
struct arg1_params{
	double Enu;	//neutrino energy
};


double pi = M_PI;
const double GF = 1E-5;	//GeV^-2
const double M = 1.0;		//isoscalar effective mass in GeV
const double Mw = 100.0;		//GeV
double q2max = pdf->q2Max();
double q2min = pdf->q2Min();

/* function that returns the differential cross-section*/
double funcv(double *x, size_t dim, void * p){
	
	/*parameters*/
	struct arg1_params * fp = (struct arg1_params*) p;
	double Enu = fp->Enu;
	double Q2 = 2 * M * Enu * x[0] * x[1];
//	double Q2 = Mw*Mw;
	double u = pdf->xfxQ2(2, x[0], Q2);
	double d = pdf->xfxQ2(1, x[0], Q2);
	double ubar = pdf->xfxQ2(-2, x[0], Q2);
	double dbar = pdf->xfxQ2(-1, x[0], Q2);
	/*distributions*/
	double xq = (u-ubar + d-dbar)/2.0;
	double xqbar = 0.0;

	/*differential cross section*/
	double eta = (2 * M * Enu/pi) * pow(GF, 2) * pow( ( Mw*Mw/(Q2 + Mw*Mw) ), 2 );
	double f = eta * ( xq + xqbar * pow( (1.0 - x[1]), 2) ) ;
	return f;
}

double funcuds(double *x, size_t dim, void * p){
	
	/*parameters*/
	struct arg1_params * fp = (struct arg1_params*) p;
	double Enu = fp->Enu;
	double Q2 = 2 * M * Enu * x[0] * x[1];
//	double Q2 = Mw*Mw;
	double ubar = pdf->xfxQ2(-2, x[0], Q2);
	double dbar = pdf->xfxQ2(-1, x[0], Q2);

	/*distributions*/
	double xq = (ubar + dbar)/2.0;
	double xqbar = (ubar + dbar)/2.0 ;

	/*differential cross section*/
	double eta = ( 2 * M*Enu/pi ) * pow(GF, 2) * pow( ( Mw*Mw/(Q2 + Mw*Mw) ), 2 );
	double f = eta * ( xq + xqbar * pow( (1.0 - x[1]), 2) ) ;
	return f;
}
double funcs(double *x, size_t dim, void * p){
	
	/*parameters*/
	struct arg1_params * fp = (struct arg1_params*) p;
	double Enu = fp->Enu;
	double Q2 = 2 * M * Enu * x[0] * x[1];

	/* Make a PDF for each quark/anti-quark type. PIDs are:
	down	up	strange	charm	beauty	truth
	1	2	3	4	5	6 */	
	double s = pdf->xfxQ2(3, x[0], Q2);
	/*distributions*/
	double xq = s ;
	double xqbar = 0.0;

	/*differential cross section*/
	double eta = ( 2 * M*Enu/pi ) * pow(GF, 2) * pow( ( Mw*Mw/(Q2 + Mw*Mw) ), 2 );
	double f = eta* ( xq + xqbar * pow( (1.0 - x[1]), 2) ) ;
	return f;
}
double funcc(double *x, size_t dim, void * p){
	
	/*parameters*/
	struct arg1_params * fp = (struct arg1_params*) p;
	double Enu = fp->Enu;
	double mparton = 1.2;
	double Q2 = 2.0 * M * Enu * x[0] * x[1];

	double cbar = pdf->xfxQ2(-4, x[0], Q2);
	/*distributions*/
	double xq = 0.0;
	double xqbar = cbar;

	/*differential cross section*/
	double eta = ( 2 * M*Enu/pi ) * pow(GF, 2) * pow( ( Mw*Mw/(Q2 + Mw*Mw) ), 2 );
	double f = eta * ( xq + xqbar * pow( (1.0 - x[1]), 2) ) ;
	if (Enu < 100.0){
		f = 0.0;
	}
	return f;
}
double funcb(double *x, size_t dim, void * p){
	
	/*parameters*/
	struct arg1_params * fp = (struct arg1_params*) p;
	double Enu = fp->Enu;
	double mparton = 173;
	double Q2 = 2 * M * Enu * x[0] * x[1] - pow(mparton, 2);
	if (Q2 < q2min){
		Q2 = q2min;
	}

	double b = pdf->xfxQ2(5, x[0], Q2);
	/*distributions*/
	double xq = 0.0;
	double xqbar = b;

	/*differential cross section*/
	double eta = ( 2 * M*Enu/pi ) * pow(GF, 2) * pow( ( Mw*Mw/(Q2 + Mw*Mw) ), 2 );
	double f = eta* ( xq + xqbar * pow( (1.0 - x[1]), 2) ) ;
	return f;
}



int main(){

	/*use these to time the program*/
	clock_t ti, tf;
	ti = clock();
	
	cout << "initializing... " << '\n';

	/* relevant variables */
	double Sigmaccv = 0.0;
	double Sigmaccuds = 0.0;
	double Sigmaccs = 0.0;
	double Sigmaccc = 0.0;
	double Sigmaccb = 0.0;
	double Err = 0.0;
	size_t calls = 100000;	//# function calls for VEGAS
	
	/* initial probe energy in GeV*/
	double Eneutrino = 10.0;

	/* integration bounds for {x, y} */
	const double xlow = 1E-8;
	const double ylow = 0.0;
	const double xup = pdf->xMax();
	const double yup = 1.0;
	double lower[2] = {xlow, ylow};	
	double upper[2] = {xup, yup};
	cout << '\n' << "Q2min: " << q2min << '\n' << "Q2max: " << q2max << '\n' << "xmin: " << xlow << '\n' << "xmax: " << xup << '\n';
	/*initialize PDF stuff */
	setVerbosity(SILENT);
	initPDFSet(setname, setmember);	

	/*set up VEGAS stuff and random number generator*/
	gsl_monte_function Fv;
	gsl_monte_function Fuds;
	gsl_monte_function Fs;
	gsl_monte_function Fc;
	gsl_monte_function Fb;
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	gsl_rng *r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	/*initialize functions for integration */
	Fv.f = &funcv;
	Fuds.f = &funcuds;
	Fs.f = &funcs;
	Fc.f = &funcc;
	Fb.f = &funcb;
	Fv.dim = 2;
	Fuds.dim = 2;
	Fs.dim = 2;
	Fc.dim = 2;
	Fb.dim = 2;

	double comparison[13];
	comparison[0] = 7.77e-38;
	comparison[1] = 6.97e-37;
	comparison[2] = 6.25e-36;
	comparison[3] = 4.54e-35;
	comparison[4] = 1.96e-34;
	comparison[5] = 6.11e-34;
	comparison[6] = 1.76e-33;
	comparison[7] = 4.78e-33;
	comparison[8] = 1.23e-32;
	comparison[9] = 3.01e-32;
	comparison[10] = 7.06e-32; 
	comparison[11] = 1.59e-31;
	comparison[12] = 0.0;
	

	/* open files to write */
	ofstream myfile1;
	ofstream myfile2;
	myfile1.open("NuN1.dat");
	myfile2.open("NuN2.dat");
	int i = 0; //dummy counter

	/* START INTEGRATING, keep cross-section as function of probe energy */
	while (Eneutrino <= 1E13){
		
		/*integration workspace*/
		gsl_monte_vegas_state * s =  gsl_monte_vegas_alloc(2);
		struct arg1_params P = {Eneutrino};
		Fv.params = &P;

		Fuds.params = &P;
		
		Fs.params = &P;

		Fc.params = &P;

		Fb.params = &P;

		/* check status */
		if (i == 0){
			cout << "working..." << '\n';
		}
		
		/*integration over Bjorken variables*/
		gsl_monte_vegas_integrate(&Fv, lower, upper, 2, calls, r, s, &Sigmaccv, &Err);
		gsl_monte_vegas_integrate(&Fuds, lower, upper, 2, calls, r, s, &Sigmaccuds, &Err);
		gsl_monte_vegas_integrate(&Fs, lower, upper, 2, calls, r, s, &Sigmaccs, &Err);
		gsl_monte_vegas_integrate(&Fc, lower, upper, 2, calls, r, s, &Sigmaccc, &Err);
		gsl_monte_vegas_integrate(&Fb, lower, upper, 2, calls, r, s, &Sigmaccb, &Err);

		/* write to files */
		myfile1 <<(Eneutrino) << " " << log10(Sigmaccv/Eneutrino *conv) << " " <<  log10(Sigmaccuds/Eneutrino *conv) << " " << log10(Sigmaccs/Eneutrino *conv)  << " " << log10(Sigmaccc/Eneutrino *conv)  << " " << log10(Sigmaccb/Eneutrino *conv) << " " << log10( (Sigmaccv + Sigmaccuds + Sigmaccs + Sigmaccc + Sigmaccb)/Eneutrino *conv) << '\n' ;

		myfile2 <<(Eneutrino) << " " << log10(Sigmaccv *conv2) << " " <<  log10(Sigmaccuds *conv2) << " " << log10(Sigmaccs*conv2)  << " " << log10(Sigmaccc *conv2)  << " " << log10(Sigmaccb *conv2) << " " << ( (Sigmaccv + Sigmaccuds + Sigmaccs + Sigmaccc + Sigmaccb) *conv2) << " " << comparison[i] << '\n' ;

		/* set probe energy, perpare to start again */
		Eneutrino = Eneutrino * 10.0;
		i++;
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


# deepinelasticscattering
