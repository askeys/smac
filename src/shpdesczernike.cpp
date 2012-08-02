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
 *  zerndesc.cpp
 *  libsmac
 *
 *  Created by askeys on 2/5/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "shpdesczernike.h"
#include "shpdescfourier.h"
#include "CellList.h"
#include "preprocess.h"
#include "space.h"
#include <sstream>
#include <iostream>


#include <assert.h>
#include <cstdlib>


namespace smac
{
	double plgndr(int l, int m, double x);

    /**
    \brief recursive factorial routine
    \param n is number to compute factorial of
    \note routine takes and returns doubles, since I was running into large 
        number problems for int.  Maybe there is a better solution...
    */
	static double factorial(double n)
	{ 
		return (n < 2) ? 1 : n*factorial(n - 1); 
	}

	/**
	\brief computes a 2d zernike moment \f$V_{nm}\f$
	\param shape is the shape data (coordinates / weights)
	\param sd is the output descriptor
	\param arg is a pointer to a zernike_info struct

	This function is not typically used directly.  It should be called by
	"zernikedesc{dim}d.
	*/
	void zernmom2d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		coordlist_t &x = shape.x;
		weight_t& f = shape.f;

		zernike_info* zarg = static_cast<zernike_info*>(arg);
		int n = zarg->_n;
		int m = zarg->_m;

		int nradterm = 1+((n-m)/2);

		std::vector<double> radcoeffs(nradterm, 0.0);
		std::vector<double> radpowers(nradterm, 0);
        

		assert( n >= 0);
		assert( m <= n );
		// Make sure n-m is not odd:
		assert( ( ( n-m) & 1 ) == 0 );

		int sign = -1;

		for ( int s = 0; s < nradterm; s++ ) {
			sign *= -1;
			radpowers[s] = n - 2 * s;
			radcoeffs[s]  = sign * factorial( n - s ) /
							(factorial( s ) * factorial( ( n + m ) / 2 - s ) *
							factorial( ( n - m ) / 2 - s ) );
		}


		component_t zmn = component_t(0.0, 0.0);

		for (unsigned int i=0; i<x.size(); i++) {

			double r = sqrt (x[i][0]*x[i][0] + x[i][1]*x[i][1]);
			if (r > 1.0) {
                std::cerr << "*** ERROR: zernmom2d; object not on unit disk\n";
                exit(1);
            }
			double Rnm = 0.0;
			for (int s = 0; s < nradterm; s++) {
				Rnm += radcoeffs[s] * pow( r, radpowers[s] );
			}

			double theta = atan2(x[i][1], x[i][0]);
            //std::cerr << r << "\t" << n << "\t" << m << "\t" << Rnm/(2*(n+1)) << std::endl;
			component_t Vnm_conj = f[i]*Rnm*component_t(cos(m*theta), -sin(m*theta));
			zmn += Vnm_conj;
		}

		sd.clear();
		zmn *= (n+1) / M_PI;
		sd.push_back(zmn/(double)x.size());
	}


	/**
	\brief computes a 3d zernike moment \f$Omega_{nl}\f$
	\param shape is the shape data (coordinates / weights)
	\param sd is the output descriptor
	\param arg is a pointer to a zernike_info struct
	*/
	void zernmom3d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		coordlist_t &x = shape.x;
		weight_t& f = shape.f;
		zernike_info* zarg = static_cast<zernike_info*>(arg);
		int n = zarg->_n;
		int l = zarg->_l;
		shpdesc_t omega_nlm;
		omega_nlm.resize(2*l+1);
		for (int i=0; i<2*l+1; i++) {
			omega_nlm[i] = component_t(0.0, 0.0);
		}

		int nradterm ( 1 + ( ( n - l ) / 2 ) );
		std::vector<double> radcoeffs(nradterm, 0.0);
		std::vector<double> radpowers(nradterm, 0);

		assert(n >= 0);
		assert(l <= n );
		// Make sure n-m is not odd:
		assert((( n-l) & 1 ) == 0 );

        //actually this starts at 1, see below
		int sign = -1;

		for ( int s = 0; s < nradterm; s++ ) {
			sign *= -1;
			radpowers[s] = n - 2 * s;
			radcoeffs[s]  = sign * factorial(n - s) / 
							(factorial(s) * factorial((n+l)/2-s)* 
							factorial((n-l)/2-s));
		}

		for (unsigned int i=0; i<x.size(); i++) {

			//calculate the radial polynomial
			double r = sqrt (
				x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2]);

			if (r > 1.0) {
                std::cerr << "*** ERROR: zernmom3d; object not on unit ball\n";
                exit(1);
            }


			double Rnl = 0.0;
			for (int s = 0; s < nradterm; s++) {
				Rnl += radcoeffs[s] * pow( r, radpowers[s] );
			}

			//calculate the spherical harmonics
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

			double c, ratio_of_two_factorial_terms = 1.0;

			std::vector< std::complex<double> > Ylm(2*l+1, 0.0);
			for (int m=0; m<=l; m++) {
				if (m > 0) {
					ratio_of_two_factorial_terms *= (l + m) * (l - m + 1);
				}
				c = sqrt((2.0*l+1.0) /(4.0*M_PI)/ratio_of_two_factorial_terms);
				c *= plgndr(l, m, cos_theta);
				component_t ylm(c * cos(m*phi), c * sin(m*phi));
				Ylm[m+l] += ylm;

				if (m > 0) {
					int sign = (int)(pow(-1.0 , (double)m));
					component_t ylminusm(
						sign*c*cos(m*phi), -1.0*sign*c*sin(m*phi)
					);
					Ylm[-m+l] += ylminusm;
				}
			}

			for (unsigned int ii=0; ii<Ylm.size(); ii++) {
				 omega_nlm[ii] += f[i]*Rnl*Ylm[ii];
			}
		}

		double inv_size = 1.0 / (double)x.size();
		for (unsigned int ii=0; ii<2*l+1; ii++) {
			omega_nlm[ii] *= 3.0 / (4.0 * M_PI) * inv_size;
		}

		if (zarg->invariant) {
			sd.push_back(component_t(invariantq(omega_nlm), 0.0));
		}
		else {
			sd = omega_nlm;
		}

	}

	/**
	\brief the 2d zernike descriptor
	\ingroup shpdesc
	\param shape is the shape data (coordinates / weights)
	\param sd is the output descriptor
	\param arg is a pointer to a zernike_info struct
	*/
	void zerndesc2d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		zernike_info* zarg = static_cast<zernike_info*>(arg);
		std::vector<int>& moment = zarg->moment;

		sd.clear();

        //normalize the shape, if necessary:
        if (zarg->normmethod == RADIALNORMCENTROID ||
            zarg->normmethod == RADIALEXCLUDECENTROID) {
            transcentroid(shape.x);
        }

        if (zarg->normmethod == RADIALNORM ||
            zarg->normmethod == RADIALNORMCENTROID) {
            assert(zarg->rmax < 1.0);
            transcentroid(shape.x);
            normuball(shape.x);
            rescale(shape.x, zarg->rmax);
        }
        else if (zarg->normmethod == RADIALEXCLUDE
            || zarg->normmethod == RADIALEXCLUDECENTROID) {
            shapedata_t shapetmp(shape.x, shape.f);
            shape.x.clear();
            shape.f.clear();
            coord_t origin(3, 0.0);
            double rsq = zarg->rmax*zarg->rmax;
            for (unsigned int i=0; i<shapetmp.x.size(); i++) {
                if (distancesq(shapetmp.x[i], origin) < rsq) {
                    shape.x.push_back(shapetmp.x[i]);
                    shape.f.push_back(shapetmp.f[i]);
                }
            }
        }


		for (unsigned int i=0; i<=moment.size(); i++) {
			int n = moment[i];
			for (int m=n%2; m<=n; m+=2) {
                if (n > 100) {
                    break;
                }
				zarg->_n = n;
				zarg->_m = m;
				shpdesc_t sd_temp;
				zernmom2d(shape, sd_temp, arg);
				sd.push_back(sd_temp[0]);
			}
		}

		if (zarg->invariant) {
			for (int i=0; i<sd.size(); i++) {
				sd[i] = component_t(
					sd[i].real()*sd[i].real() + sd[i].imag()*sd[i].imag(), 0.0);
			}
		}
	}

	/**
	\brief the 3d zernike descriptor
	\ingroup shpdesc
	\param shape is the shape data (coordinates / weights)
	\param sd is the output descriptor
	\param arg is a pointer to a zernike_info struct
	*/
	void zerndesc3d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		zernike_info* zarg = static_cast<zernike_info*>(arg);
		std::vector<int>& moment = zarg->moment;

        sd.clear();

        //normalize the shape, if necessary:
        if (zarg->normmethod == RADIALNORMCENTROID ||
            zarg->normmethod == RADIALEXCLUDECENTROID) {
            transcentroid(shape.x);
        }

        if (zarg->normmethod == RADIALNORM ||
            zarg->normmethod == RADIALNORMCENTROID) {
            assert(zarg->rmax < 1.0);
            transcentroid(shape.x);
            normuball(shape.x);
            rescale(shape.x, zarg->rmax);
        }
        else if (zarg->normmethod == RADIALEXCLUDE
            || zarg->normmethod == RADIALEXCLUDECENTROID) {
            shapedata_t shapetmp(shape.x, shape.f);
            shape.x.clear();
            shape.f.clear();
            coord_t origin(3, 0.0);
            double rsq = zarg->rmax*zarg->rmax;
            for (unsigned int i=0; i<shapetmp.x.size(); i++) {
                if (distancesq(shapetmp.x[i], origin) < rsq) {
                    shape.x.push_back(shapetmp.x[i]);
                    shape.f.push_back(shapetmp.f[i]);
                }
            }
        }

        //compute the moments
		for (unsigned int i=0; i<moment.size(); i++) {
			int n = moment[i];
			for (int l=n%2; l<=n; l+=2) {
				zarg->_n = n;
				zarg->_l = l;
				shpdesc_t moment_nl;
				zernmom3d(shape, moment_nl, arg);
				sd.insert(sd.begin(), moment_nl.begin(), moment_nl.end());
			}
		}
	}

/*
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
*/

}
