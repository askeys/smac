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
 *  shape_distribution.cpp
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "shpdescdist.h"
#include "HistogramConst.hpp"
#include "space.h"
#include <assert.h>

namespace smac 
{
	/**
	\param x is a 2d array of coordinates.  Ordering is {double* ptr2x1, ... 
		double* ptr2xn}.  The double* ptr2xi must point to a block of memory 
		of size > 3*sizeof(double).
	\param binwidth is the bin resolution (i.e., the width of each bin in the 
		histogram)
	\param sd is an array of real numbers representing the computed shape 
		descriptor 
	*/
	void shpdistD2(shapedata_t& s, shpdesc_t& sd, arg_t arg) 
	{
		shpdistD2_info* s_arg = static_cast<shpdistD2_info*>(arg);			
        shpdistD2_info default_info;
		if (arg == 0x0) {
            s_arg = &default_info;
        }
        
        box_t& box = s_arg->box;
        coordlist_t& x = s.x;
		int n = x.size();
		double rmax = s_arg->rmax;
		double rmin = s_arg->rmin;
		double rstep = (rmax - rmin) / (double)(s_arg->nbin);
		HistogramConst<double> h(rmin, rstep, rmax);
		h.insert(0.0, 0.0);														//make sure the histogram starts off at 0.0
		double ntotal = 0.0;
		
		for( int i=0; i<n; i++) {
			for( int j=i+1; j<n; j++) {
                double r = distance(x[i], x[j], box.period, box.periodic);
                if (r >= rmin && r < rmax) {
                    h.insert(r, 2.0);
                    ntotal += 2.0;
                }
			}
		}
		
		sd.clear();
		int size = h.getNumberOfBins();
		for (int i=0; i<size; i++) {
			std::pair<double, double> bini = h.getBin(i);
			sd.push_back(component_t(bini.second / ntotal));
		}		
	}

	/**
	\param x is a 2d array of coordinates.  Ordering is {double* ptr2x1, ... 
		double* ptr2xn}.  The double* ptr2xi must point to a block of memory 
		of size > 3*sizeof(double).
	\param binwidth is the bin resolution (i.e., the width of each bin in the 
		histogram)
	\param sd is an array of real numbers representing the computed shape 
		descriptor 
	*/
	void shpdistA3(shapedata_t& s, shpdesc_t& sd, arg_t arg)
	{
		shpdistA3_info* s_arg = static_cast<shpdistA3_info*>(arg);			
        shpdistA3_info default_info;
		if (arg == 0x0) {
            s_arg = &default_info;
        }
        box_t& box = s_arg->box;
        coordlist_t& x = s.x;
		int n = x.size();
		int nbin = s_arg->nbin;
		double binwidth = M_PI / (double)nbin;
		
		HistogramConst<double> h(0.0, binwidth, nbin*binwidth);	
		double ntotal = 0.0;
		
		//algorithm double counts, but properly includes all angles
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {
				if (i != j) {
					double dx = x[i][0] - x[j][0]; 
                    pbc(dx, box.period[0], box.periodic[0]);
					double dy = x[i][1] - x[j][1]; 
                    pbc(dy, box.period[1], box.periodic[1]);
					double dz = x[i][2] - x[j][2];
                    pbc(dz, box.period[2], box.periodic[2]);
					double rij = dx*dx + dy*dy + dz*dz;

					for (int k=0; k<n; k++) {
						if (i!=k  && j!=k) {
							double dx = x[j][0] - x[k][0]; 
                            pbc(dx, box.period[0], box.periodic[0]);
							double dy = x[j][1] - x[k][1]; 
                            pbc(dy, box.period[1], box.periodic[1]);
							double dz = x[j][2] - x[k][2];
                            pbc(dz, box.period[2], box.periodic[2]);
							double rik = dx*dx + dy*dy + dz*dz;
							
							double x1 = x[i][0] - x[j][0];
                            pbc(x1, box.period[0], box.periodic[0]);
							double y1 = x[i][1] - x[j][1];
                            pbc(y1, box.period[1], box.periodic[1]);
							double z1 = x[i][2] - x[j][2];
                            pbc(z1, box.period[2], box.periodic[2]);

							double x2 = x[k][0] - x[j][0];
                            pbc(x2, box.period[0], box.periodic[0]);
							double y2 = x[k][1] - x[j][1];
                            pbc(y2, box.period[1], box.periodic[1]);
							double z2 = x[k][2] - x[j][2];
                            pbc(z2, box.period[2], box.periodic[2]);

							double dot = x1*x2 + y1*y2 + z1*z2;					
							double theta = acos(dot/sqrt(rij*rik));
							h.insert(theta);
							ntotal += 1.0;
						}
					}
				}
			}
		}
		/*need to divide by half since we double count*/
		sd.clear();
		int size = h.getNumberOfBins();
		for (int i=0; i<size; i++) {
			std::pair<double, double> bini = h.getBin(i);
			sd.push_back(component_t(bini.second/ntotal));
		}		
	}
	
	/**
	\brief computes the probabiltiy distribution function for an arbitrary 
		shape descriptor
	\param sdlist is a list of shape descriptors from which to compute the 
		probabilty distribution function
	\param sd is the output shape descriptor
	\param arg is a pointer to list of arguments of type shphist_info
	*/
	void shpdist(shpdesclist_t& sdlist, shpdesc_t& sd, arg_t arg)
	{
		shpdist_info* s_arg = static_cast<shpdist_info*>(arg);			
		assert(s_arg != NULL);														

		double binmax = s_arg->max;
		double binmin = s_arg->min;
		double binstep = (binmax - binmin) / (double)(s_arg->nbin);
		        
		assert(sdlist.size() > 0);
		int sdsize = sdlist[0].size();
		sd.clear();
				
		for (int j=0; j<sdsize; j++) {
			HistogramConst<double> h_real(binmin, binstep, binmax);	
			HistogramConst<double> h_imag(binmin, binstep, binmax);
			int ntotal = 0;	
			for (unsigned int i=0; i<sdlist.size(); i++) {
				h_real.insert(sdlist[i][j].real());
				h_imag.insert(sdlist[i][j].imag());
				ntotal++;
			}
            
            //std::cerr << binmin << "\t" << binmax << "\t" << binstep << "\t" << h_real.getNumberOfBins() << std::endl;
			int nbins = h_real.getNumberOfBins();
			for (int k=0; k<nbins; k++) {
                std::pair<double, double> binreal = h_real.getBin(k);
				std::pair<double, double> binimag = h_imag.getBin(k);
				if (s_arg->comptype == REAL) {
                    //std::cerr<<k << "\t" << binreal.second/(double)ntotal << std::endl;
                    sd.push_back(
                        component_t(binreal.second/(double)ntotal, 0.0));
                }
                else {
                    sd.push_back(component_t(binreal.second/(double)ntotal, 
                        binimag.second/(double)ntotal));
                }
			}		
		}
	}
}
