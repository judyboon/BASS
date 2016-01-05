/*
  Bayesian group factor Analysis with Structured Sparsity (BASS)
  Copyright (C) 2015 Shiwen Zhao
  
  This source code is downloaded from 
  http://jonaswallin.github.io/articles/2013/07/simulation-of-gig-distribution/
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "gig_par.h"
#include <iostream>

// should be made parallell?
//function [Eig,Eg, Elogg] = Egig(p, a, b)
void Egig(const int n,const double* p, const double* a, const double* b, double* Eg , double* Eig , double* Elogg )
{
	
	int i;
	double p_, p1_, sqrt_ab;
	double log_Kp, log_Kp1;
	
	for( i = 0; i < n; i++)
	{
		p_ = fabs(p[i]);
		
		sqrt_ab = sqrt(a[i] * b[i]);
		log_Kp  = gsl_sf_bessel_lnKnu(p_ , sqrt_ab);
		if(Eg){
			p1_ = fabs( p[i] + 1);
			log_Kp1 = gsl_sf_bessel_lnKnu(p1_, sqrt_ab);
			Eg[i] = exp( (log(b[i]) - log(a[i])) / 2 + log_Kp1 -log_Kp);
		}
		
		if(Eig)
		{
			p1_ = fabs( p[i] - 1);
			log_Kp1 = gsl_sf_bessel_lnKnu(p1_, sqrt_ab);
			Eig[i] = exp( ( - log(b[i]) + log(a[i])) / 2 + log_Kp1 -log_Kp);
		}
		if(Elogg){
			double eps = 1e-7;//fabs(p[i] ) *1.4901e-08;
			p1_ = fabs( p[i] + eps);
			log_Kp1 = gsl_sf_bessel_lnKnu(p1_, sqrt_ab);
			Elogg[i] =  (log(b[i] ) -log(a[i])) /2; 
			Elogg[i] += (exp(log_Kp1 -log_Kp  ) - 1) / eps;
		}
	}
} 

double gig::sqrt_gig_ratio( double x_in, double m_in){
	
	double ret = pow(x_in / m_in, (p_abs - 1)/2);
	ret *= exp( ( m_in + 1/m_in - x_in - 1/x_in  )/(2 * two_d_beta ) );
	return ret;
}	

void gig::algortihm_region1()
{
	two_d_beta = 2. / beta;
	
	// calc mode
	mode = sqrt(pow(p_abs - 1,2) + pow(beta, 2) ) + (p_abs - 1);
	mode /= beta;
	
	
	// COMPUTING minimal bounding rectangle through
	//Setting up and solving Caradanos to get x_m,x_p
	c1 = - two_d_beta * (p_abs + 1)  - mode;
	c2 =   two_d_beta * (p_abs - 1)  * mode - 1;
	c3 = c2 - pow(c1, 2) / 3.;
	c4 = 2 * pow(c1, 3)/ 27. - c1 * c2 / 3 + mode;
	c5 = acos(- c4/2 * sqrt(- 27/ pow( c3, 3 )));
 	c6 = sqrt(-4 * c3 / 3. );
	x_m = c6 * cos( c5 / 3 + 4 * M_PI / 3 ) - c1 / 3;
	x_p = c6 * cos(c5 / 3               ) - c1 / 3;
	
	
	// triangle to simulate
	u_m_div_v = (x_m - mode) * sqrt_gig_ratio(x_m, mode);
	u_p_div_v = (x_p - mode) * sqrt_gig_ratio(x_p, mode);
	u_p_m = u_p_div_v - u_m_div_v;
	
	// generating
	bool new_value = true;
	do{
		U = RngSt->RandU01();
		V = RngSt->RandU01();
		x  = u_p_m * (U / V) + (u_m_div_v / V)  + mode; 
		if(x < 0)
			continue;
		
		if(V <= sqrt_gig_ratio(x, mode))
			new_value = false;
			
	}while(new_value );
}

void gig::algortihm_region2()
{
	two_d_beta = 2. / beta;
	// calc mode
	mode =  beta;
	mode /= sqrt(pow(1 - p_abs ,2) + pow(beta, 2) ) + (1 - p_abs );	
	//std::cout << "mode =" << mode <<"\n";
	// benerating constants and domains
	x0 = beta / (1 -p_abs);
	xs = std::max(x0, two_d_beta);
	
	c1 = gig_propto(mode);
	A1 = c1 * x0;
	if( x0 < two_d_beta)
	{
		c2 = exp(- beta);
		if(p_abs >0){
			x0_pow_p = pow(x0, p_abs); 
			A2 = c2 *( pow(two_d_beta, p_abs) - x0_pow_p)/ p_abs;
		}
		else // |p|=0
			A2 = c2 * log(two_d_beta / beta);
	}else{
		c2 = 0;
		A2 = 0;
	}
	c3 = pow(xs, p_abs - 1);
	A3 = 2 * c3 * exp(-xs / two_d_beta) / beta;
	A = A1 + A2 + A3;
	// generator
	do{
		U =  RngSt->RandU01();
		V =  RngSt->RandU01() * A;
		
		if(V <= A1){ // region (0, x0)
			x = x0 * V / A1;
			c4 = c1;
		}else if(V <= A1 + A2){ // region (x0, two_d_beta)
			V = V - A1;
			if(p_abs > 0)
				x = pow(x0_pow_p + V * p_abs / c2, 1. / p_abs);
			else
				x = beta * exp(V * exp(beta));
			c4 = c2 * pow(x, p_abs - 1);	
		}else{ // region (two_d_beta, \infinty)
			V = V - A1 - A2;
			x = - two_d_beta * log( exp( - xs / two_d_beta) - V  / (c3 *two_d_beta ));
			c4 = c3 * exp( - x / two_d_beta);
		}
	}while(U * c4 >= gig_propto(x));
}
 
void gig::algortihm_region3()
{
	two_d_beta = 2. / beta;
	// calc mode
	mode =  beta;
	mode /= sqrt(pow( 1 - p_abs ,2) + pow(beta, 2) ) + 1 - p_abs ;	
	// computing miniminal bounding rectangle
	x_p = 1 + p_abs + sqrt( pow(1 + p_abs,2) + pow(beta, 2) );
	x_p /= beta;
	v_p = sqrt_gig_propto(mode);
	u_p = x_p * sqrt_gig_propto(x_p);
	int count = 0;
	do{
		U =  RngSt->RandU01() * u_p;
		V =  RngSt->RandU01() * v_p;
		x = U / V;
	}while(V >= sqrt_gig_propto(x));
}

void gig::rgig(double& x_out, const double  *p,const double *a, const double *b) {
	
	
	// comment if b = 0 should impliment different algorithm

	 alpha = sqrt(*a / *b);
	 beta  = sqrt(*a * *b);
	// reformulating the model to:
	//  \propto alpha^{-p} x^{p-1} e^{-beta/2 * ( (x/alpha)^-1 + (x/alpha)) }
	
	// utilize that x^-1 \propto GIG(-p,\alpha,\beta) 
	p_abs = fabs(*p);
	if( beta > 1 | p_abs > 1){
		//std::cout << "reagion 1\n";
		algortihm_region1();
	}else if( beta <= std::min(1/2.,2./3. * sqrt(1-p_abs))){
		//std::cout << "reagion 2\n";
		algortihm_region2();
	}else{
		//std::cout << "reagion 3\n";
		algortihm_region3();
	}
	
	if(*p < 0)
		x_out = 1. /x;
	else
		x_out  = x;
	x_out /= alpha;
}


	
