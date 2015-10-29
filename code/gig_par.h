#ifndef GIG_PAR_H
#define GIG_PAR_H
// random generalized inverse Gaussian generator
// lambda, chi, psi using wiki notation
// *double where to put the r.v
// length of random variableser
// adding seed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "RngStream.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_poly.h>


#define Calloc(n,type) (type *)calloc((size_t)(n),sizeof(type))



/*
 * Computes the Expcation of GIG ,GIG^-1 and log(GIG)
 *  @param n number of observations
	@param p parameter vector for GIG  
	@param a parameter vector for GIG
	@param b parameter vector for GIG
	@param Eg vector where storing E[GIG], (if 0 dont calculate)
	@param Eig vector where storing E[GIG^-1], (if 0 dont calculate)
	@param Elogg vector where storing E[log(GIG)], (if 0 dont calculate)
 * 
 */ 
void Egig(const int n, const double* p, const double* a, const double* b, double* Eg=0, double* Eig = 0, double* Elogg = 0);

/*
 * Class for generating GIG random variable
 * implimented for openmp
 * 
 * pdf
 *  f(x) = \frac{(a/b)^{\lambda/2}}{2K_{\lambda}(\sqrt{ab})} x^{\lambda-1} e^{-(ax + bx^{-1})/2}
 * 
 * Allthough GIG is defined for b = 0, the algorithm is not implimented for that case
 * The algorithm comes from:
    Generating generalized Inverse Gaussian Random Variates, Hörmann and Leydold, 2013
	
	The rejction constant is unfiormly bounded from above with 2.72604 * 
	 * 
 */ 
class gig
{
private:


	
	
	double mode; // parameters for model 
	double c1, c2, c3, c4, c5, c6;
	double A, A1, A2, A3;
	double x_m, x_p;
	double v_p, u_m, u_p, u_p_m;
	double u_m_div_v, u_p_div_v;
	double x0, xs;
	double two_d_beta, x0_pow_p;
	double U, V;
	
	
	double alpha, beta, p_abs; //parameter modification
	double x;
	RngStream* RngSt; // U[0,1] simulator
public:
	gig(){ RngSt = 0; };
	~gig(){};

	/*
	 * Sett Random stream
	 */ 
	void setRngStream( RngStream* in) { RngSt = in;}

	
	
	/*
	 * sampling gig random variable
	 * @param the value where to store x
	 * @param p
	 * @param a
	 * @param b
	 */
	void rgig(double& ,const double* ,const  double* ,const double* );
	
	/*
	 * Computes x^{pa - 1} exp( - \beta/2 (x + 1/x) )
	 */ 
	double gig_propto( double x_in ){ return exp( (p_abs-1) * log(x_in)- (x_in + 1 / (x_in) )/two_d_beta); }
	
		/*
	 * Computes x^{(pa - 1)/2} exp( - \beta/4 (x + 1/x) )
	 */ 
	double sqrt_gig_propto( double x_in ){ return exp((p_abs-1)/2 *log(x_in) - (x_in + 1 / (x_in) )/(2 * two_d_beta ) ); }
		
		
	/* Computes the ratio between to square root of two gig distribtuion
	 * @ param x_in the denominator term 
	 * @ param m_in the  numerator
	 * Computes sqrt_gig_propto(x_in) / sqrt_gig_propto(m_in)
	 */ 
	double sqrt_gig_ratio( double , double );	
	
	
	
	/*
		Algorithm if the parameter region is
		sqrt(a * b)> 1 or abs(p) > 1
		Algorithm 3 from Hörmann Leydold (Created by Dagpunar)
		
		@out stores the sampled value in GIG.x
	*/ 
	void algortihm_region1();
	/*
		Algorithm if the parameter region is
		sqrt(a * b)< = 2/3 * \sqrt{1 - \lambda} or abs(p) a <= 1
		Algorithm 1 from Hörmann Leydold 
		
		@out stores the sampled value in GIG.x
	*/ 
	void algortihm_region2();
	/*
		Algorithm if the parameter region is
		1 >) sqrt(a * b) >  2/3 * \sqrt{1 - \lambda} or abs(p) a <= 1
		Algorithm 2 from Hörmann Leydold 
		
		@out stores the sampled value in GIG.x
	*/ 
	void algortihm_region3();	
};
#endif /* RMATH_H */
