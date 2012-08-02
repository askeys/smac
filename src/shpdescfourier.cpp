/**
SMAC -- Shape Matching Analysis Code
(SMAC) Open Source Software License Copyright 2012 The Regents of the 
University of Michigan All rights reserved.

SMAC may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

You may redistribute, use, and create derivate works of SMAC, in source
and binary forms, provided you abide by the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer both in the code and
prominently in any materials provided with the distribution.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* Apart from the above required attributions, neither the name of the copyright
holder nor the names of LibTPS's contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*//*
 *  fourier_descriptor.cpp
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "shpdescfourier.h"
#include "shpdeschist.h"
#include "shpdescwigner3j.h"
#include "io.h"
#include <iostream>
#include <assert.h>


namespace smac
{
	double plgndr(int l, int m, double x);

#pragma mark FOURIER COEFFICIENTS    
	/**
	\brief computes the 2d fourier coefficient $f c_\ell \$f
	\param shape is shape data (coordinates / weights)
	\param arg is a pointer to an argument struct of type fourier_info   
	\param sd is an array of complex numbers representing the computed shape 
		descriptor 
	*/
	void fouriercoeff2d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		coordlist_t& x = shape.x;
		weight_t& f = shape.f;
		fourier_info* f_arg = static_cast<fourier_info*>(arg);			
		assert(f_arg != NULL);														
		int k = f_arg->_frequency;
		std::vector<double> shells = f_arg->shells;
		int nshells = shells.size()-1;
		assert(nshells >= 1);
		assert(k >= 0);

        if (f_arg->trigtablesize > 0) {
            int table_size = f_arg->trigtablesize;
            if (f_arg->_costable.size() != (unsigned int)table_size) {
                f_arg->_costable.resize(table_size);
                f_arg->_sintable.resize(table_size);
                for (int i=0; i<table_size; i++) {
                    f_arg->_costable[i] = cos((double)i/(double)(table_size-1)*2.0*M_PI);
                    f_arg->_sintable[i] = sin((double)i/(double)(table_size-1)*2.0*M_PI);
                    //std::cerr << i << "\t" << f_arg->_costable[i] << "\t" << f_arg->_sintable[i] << "\n";
                }
            }
        }
        
		std::vector<double> rsq = f_arg->shells;
		for (unsigned int i=0; i<rsq.size(); i++) {
			rsq[i] *= rsq[i];
		}
        
		sd.clear();
		int n = x.size();
		for (int s = 0; s<nshells; s++) {
			int ns=0;
			component_t c(0.0, 0.0);
			for (int i=0; i<n; i++) {
				double rsqi = x[i][0]*x[i][0] + x[i][1]*x[i][1];
				if (rsqi >= rsq[s] && rsqi < rsq[s+1]) {
					double theta = atan2(x[i][1],x[i][0]) + M_PI;				
                    c += f[i]*component_t(cos(k*theta), -sin(k*theta));
					ns++;
				}
			}
			
			if (ns > 0) {
				c /= (double)ns;
			}            
			if (f_arg->invariant) {
				c = abs(c);
			}
			sd.push_back(c);
		}
	}

		
	/**
	\brief computes the 3d fourier coefficient $f q_\ell \$f
	\param shape is shape data (coordinates / weights)
	\param arg is a pointer to an argument struct of tye fourier_info   
	\param sd is an array of complex numbers representing the computed shape 
		descriptor 
	*/
	void fouriercoeff3d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		coordlist_t& x = shape.x;
		weight_t& f = shape.f;
		fourier_info* f_arg = static_cast<fourier_info*>(arg);			
		assert(f_arg != NULL);														
		int l = f_arg->_frequency;
		std::vector<double> r_shell = f_arg->shells;
		int nshells = r_shell.size()-1;
		sd.clear();

		int n = x.size();
		assert(l>=0);
		
		for (int s=0; s<nshells; s++) {
		
			shpdesc_t ylm_shell(2*l+1, 0.0);			
			int ns = 0;
			
			for (int i=0; i<n; i++) {			
				
			//-- STEP 1: Convert to spherical coordinates --
				double r = sqrt(x[i][0]*x[i][0] + 
								x[i][1]*x[i][1] + 
								x[i][2]*x[i][2]);
				
				if (r >= r_shell[s] && r < r_shell[s+1]) {
					
					ns++;
					double theta, phi;
					
					if(x[i][0] ==0 ) {
						if(x[i][1]==0) { 
							phi = 0.0;
						}
						else if(x[i][1] > 0 ) {
							phi = M_PI * 0.5;
						}
						else {  
							phi = 3.0 * M_PI * 0.5;
						}
					}
					else {
						phi = atan( x[i][1]/x[i][0] );
						if( x[i][0] < 0.0 ) {
							phi += M_PI;
						}
						else if( x[i][1] < 0.0 ) {
							phi += 2.0 * M_PI;
						}
					}

					double cos_theta = x[i][2]/r;
					theta = acos(cos_theta);
					
					
					//-- STEP 2: compute Ylm(theta, phi) for each theta, phi --
					//Notes:
					//1) We compute Yl(m) and Yl(-m) concurrently within the same loop 
					//   by using the symmetry of the ls
					//2) We don't really need to calculate factorials and then divide 
					//   for the the Ylm function, since we really only need the ratio 
					//   of factorials.  Thus, we use the ratio_of_two_factorial_terms 
					//   temp variable which a) allows us again to move the calculation 
					//   into a single loop for efficiency and b) solves the problem of 
					//   roundoff error that occurs if we divide two huge factorials 
					//   for large l
					
					double c, ratio_of_two_factorial_terms = 1.0;
			
					for (int m=0; m<=l; m++) {
						if (m > 0) {
							ratio_of_two_factorial_terms *= (l + m) * (l - m + 1);
						}
						c = sqrt((2.0*l+1.0) /(4.0*M_PI)/ratio_of_two_factorial_terms); 
						c *= plgndr(l, m, cos_theta);
						component_t ylm(c * cos(m*phi), c * sin(m*phi)); 
						ylm_shell[m+l] += f[i]*ylm;
										
						if (m > 0) {
							int sign = (int)(pow(-1.0 , (double)m));
							component_t ylminusm(
								sign*c*cos(m*phi), 
								-1.0*sign*c*sin(m*phi)); 
							ylm_shell[-m+l] += f[i]*ylminusm;
						}
					}
				}
			}
			if (ns != 0) {
				for(unsigned int i=0; i<ylm_shell.size(); i++) {
					ylm_shell[i] /= (double)ns;
				}
			}
			
			if (f_arg->invariant) {
				if (f_arg->invariant & INVARIANT_Q) {
					double q = invariantq(ylm_shell);
					sd.push_back(q);
				}
				if (f_arg->invariant & INVARIANT_W) {
					double w = invariantw(ylm_shell);
					sd.push_back(w);
				}
			}
			else {
				sd.insert(sd.end(), ylm_shell.begin(), ylm_shell.end());				
			}
		}
	}

#pragma mark FOURIER_DESCRIPTOR
	/**
	\brief computes the 2d fourier descriptor
	\ingroup shpdesc
	\param s is shape data (coordinates / weights)
	\param arg is a pointer to an argument struct of tye fourier_info   
	\param sd is an array of complex numbers representing the computed shape 
		descriptor 
	*/
	void fourierdesc2d(shapedata_t &s, shpdesc_t& sd, arg_t arg)
	{
		fourier_info* f_arg = static_cast<fourier_info*>(arg);			
		assert(f_arg != NULL);
		
		sd.clear();

		for (unsigned int i=0; i<f_arg->frequency.size(); i++) {
			shpdesc_t f_k;
			f_arg->_frequency = f_arg->frequency[i];
			fouriercoeff2d(s, f_k, arg);
			if (f_arg->invariant) {
				for (unsigned int i=0; i<f_k.size(); i++) {
					f_k[i] = abs(f_k[i]);
				}
			}
			sd.insert(sd.end(), f_k.begin(), f_k.end());
		}	
	}

	
	/**
	\brief computes the 3d fourier descriptor
	\ingroup shpdesc
	\param s is shape data (coordinates / weights)
	\param arg is a pointer to an argument struct of tye fourier_info   
	\param sd is an array of complex numbers representing the computed shape 
		descriptor 
	*/
	void fourierdesc3d(shapedata_t& s, shpdesc_t& sd, arg_t arg)
	{
		fourier_info* f_arg = static_cast<fourier_info*>(arg);			
		assert(f_arg != NULL);	
		
		sd.clear();
		for (unsigned int i=0; i<f_arg->frequency.size(); i++) {
			shpdesc_t y_lm;
			f_arg->_frequency = f_arg->frequency[i];
			fouriercoeff3d(s, y_lm, arg);
			sd.insert(sd.end(), y_lm.begin(), y_lm.end());
		}			
	}



#pragma mark FOURIER FUNCTIONS
	
	/**
	\brief computes the "q" invariant for a 3d fourier coefficient
	\param sd is a 3d fourier coefficent for which to compute q
	\return q for the input coefficient
	*/	
	double invariantq(shpdesc_t&sd)
	{
		//|qlm|^2 = the projection of the harmonic descriptor onto itelsef
		double qlm_sq = 0.0;
		for (unsigned int i=0; i<sd.size(); i++) {
			qlm_sq += sd[i].real()*sd[i].real() + sd[i].imag()*sd[i].imag();
		}
							 
		//normalize by 4pi/(2l+1) and take sqrt
		//the size of the sd vector is 2*l + 1;
		return sqrt(4.0*M_PI / (sd.size()) * qlm_sq );
	}
	
	/**
	\brief computes the "w3j" invariant for a 3d fourier coefficient
	\param sd is a 3d fourier coefficent for which to compute q
	\return q for the input coefficient
	*/	
	double invariantw(shpdesc_t&sd)
	{
		return wigner3jcoeff(sd);
	}

#pragma mark HELPER FUNCTIONS
	
#ifndef DOXYGEN_SHOULD_SKIP_THIS

	double plgndr(int l, int m, double x)
	{
		double fact,pll=0.0,pmm,pmmp1,somx2;
		int i,ll;

		if (m < 0 || m > l || fabs(x) > 1.0) {
			std::cerr << "Bad arguments in routine PLGNDR" << std::endl;
		}
		pmm=1.0;
		if (m > 0) {
			somx2=sqrt((1.0-x)*(1.0+x));
			fact=1.0;
			for (i=1;i<=m;i++) {
				pmm *= -fact*somx2;
				fact += 2.0;
			}
		}
		if (l == m)
			return pmm;
		else {
			pmmp1=x*(2*m+1)*pmm;
			if (l == (m+1))
				return pmmp1;
			else {
				for (ll=(m+2);ll<=l;ll++) {
					pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
					pmm=pmmp1;
					pmmp1=pll;
				}
				return pll;
			}
		}
	}
	
	#endif

}

