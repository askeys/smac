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
 *  shape_context.cpp
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "shpdesccontext.h"
#include <iostream>
#include <map>
#include "io.h"


#include <assert.h>

#ifdef _USE_HUNGARIAN
#include "../3rdparty/hungarian/hungarian.h"
#endif

namespace smac 
{
#pragma mark SHAPE CONTEXTS
	
	/**
	\param sc1 is a shape context descriptor
	\param sc2 is a shape context descriptor
	\param n is the number of points used to construct each shape context
	\param nradial is the number of radial bins
	\param nbin is the number of angular bins
	\bool verbose indicates whether to print messages to the screen
	*/
	void shape_context_reorder(
		shpdesclist_t& sc1, shpdesclist_t& sc2, 
		matchfun_t util_function,  
		int nradial, bool verbose)
	{
#ifdef _USE_HUNGARIAN
		int m = sc1.size()/nradial;
		int n = sc2.size()/nradial;
				
		std::vector<int> utility(m*n);
		
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				double total = 0.0;
				for (int k=0; k<nradial; k++) {
					total += util_function(
						sc1[i*nradial + k], sc2[j*nradial + k], 0.0, 100.0);
				}
				utility[i*m+j] = total;
			}
		}
		
		using std::map;
		
		//Exclude outliers:
		if (m > n) {
			map<int, int> map_util_to_index;
			for (int i=0; i<m; i++) {
				int total_m = 0;
				for (int j=0; j<n; j++) {
					total_m += utility[i*m+j];
				}
				map_util_to_index[total_m] = i;
			}
			int diff = m-n;
			map<int, int>::iterator z = map_util_to_index.begin();
			std::vector<bool> skip(m*nradial, false);
			for (int k=0; k<diff; k++) {
				skip[z->second] = true;
				++z;
			}
			shpdesclist_t sd_tmp_array;
			for (int i=0; i<m; i++) {
				if (!skip[i]) {
					for (int k=0; k<nradial; k++) {
						sd_tmp_array.push_back(sc1[i*nradial + k]);
					}
				}
			}
			sc1 = sd_tmp_array;
			shape_context_reorder(sc1, sc2, util_function, nradial, verbose);
			return;
		}
		else if (n > m) {
			map<int, int> map_util_to_index;
			for (int j=0; j<n; j++) {
				int total_n = 0;
				for (int i=0; i<m; i++) {
					total_n += utility[i*m+j];
				}
				map_util_to_index[total_n] = j;
			}
			int diff = n-m;
			map<int, int>::iterator z = map_util_to_index.begin();
			std::vector<bool> skip(n*nradial, false);
			for (int k=0; k<diff; k++) {
				skip[z->second] = true;
				std::cerr << z->second << std::endl;
				++z;
			}
			shpdesclist_t sd_tmp_array;
			for (int i=0; i<n; i++) {
				if (!skip[i]) {
					for (int k=0; k<nradial; k++) {
						sd_tmp_array.push_back(sc2[i*nradial + k]);
					}
				}
			}
			sc2 = sd_tmp_array;
			shape_context_reorder(sc1, sc2, util_function, nradial, verbose);
			return;
		}
		//Done excluding outliers
		
		hungarian_t prob;
		hungarian_init(&prob, &utility[0], m, n, HUNGARIAN_MAX);
		if (verbose) {
			hungarian_print_rating(&prob);
		}
		//int possible = hungarian_check_feasibility(&prob);
		//if (!possible) {
		//	std::cerr << "Assignment not possible\n";
		//	return;
		//}
		
		hungarian_solve(&prob);
		
		for (int i=0; i<m; i++) {
			for (int j=i+1; j<n; j++) {				
				if (j==prob.a[i] && i!=j) {
					for (int k=0; k<nradial; k++) {
						shpdesc_t temp = sc2[i*nradial + k];
						sc2[i*nradial+k] = sc2[j*nradial+k];
						sc2[j*nradial+k] = temp;
					}
				}
			}
		}
				
		if (verbose) {
			hungarian_print_assignment(&prob);
		}
		hungarian_fini(&prob);
#else
        int USE_HUNGARIAN = 0;
        assert(USE_HUNGARIAN);
#endif
	}
			
	/**
	\param x is an array of coordinates
	\param nbin is the number of angular bins
	\param rmin 
	*/
	void shape_context(shapedata_t& shape, shpdesclist_t& sd, 
		double rmin, double rstep, double rmax, double p, int& nradial, 
		shpdescfun_t shape_function, arg_t args)
		
	{
		assert(shape_function != NULL);
		assert(args != NULL);
/*		
		int n = x.size();
		assert(n > 0);
		int dim = x[0].size();
		
		int c = 1;																
		double step = 0;														
		std::vector<double> radius;												
		for (double r = rmin; r<=rmax; r=rmin+step) {							
			step = pow(c++, p) * rstep;											
			radius.push_back(r);												
		}																		
		radius.push_back(1e20);													
																				
		nradial = radius.size()-1;												
		assert(nradial > 0);													
																				
		for (int i=0; i<n; i++) {												
			for (int b = 0; b<nradial; b++) {									
				coordlist_t x_subset;											
				for (int j=0; j<n; j++) {
					coord_t dx(dim);
					double rsq = 0;
					for (int k=0; k<dim; k++) {
						dx[k] = x[i][k] - x[j][k];								
						rsq += dx[k]*dx[k];
					}
					double rij = sqrt(rsq);
					if(rij >= radius[b] && rij < radius[b+1]) {					
						x_subset.push_back(dx);
					}															
				}																

				shpdesc_t sd_subset;
				shape_function(x_subset, sd_subset, args);
				sd.push_back(sd_subset);				
			}																	
		}
		*/										
	}
}
