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
 *  shphist.cpp
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "shpdeschist.h"
#include <iostream>
#include "matchfun.h"

#ifdef _USE_GEOMETRY
#include "../3rdparty/geometry/geometry.H"
#endif



#include <assert.h>
#include <cstdlib>

namespace smac {


	/**
	@param x is a 2d array of coordinates.  Ordering is {double* ptr2x1, ...
		double* ptr2xn}.  The double* ptr2xi must point to a block of memory
		of size > 3*sizeof(double).
	@param nbin in the number of bins per one revolution around the circle
	@param sd is an array of real numbers representing the computed shape
		descriptor
	*/

	void shphist(shapedata_t&s, shpdesc_t& sd, arg_t arg)
	{
		sd.clear();
		if (s.x.size() == 0) {
			return;
		}
		if (s.x[0].size() == 3) {
			shphist3d(s, sd, arg);
		}
		else {
			shphist2d(s, sd, arg);
		}
	}

	void shphist2d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		shphist_info* s_arg;

        shphist_info default_info;
        if (arg == 0x0) {
            s_arg = &default_info;
        }
        else {
            s_arg = static_cast<shphist_info*>(arg);
        }

        int nsectors = s_arg->nsectors;
		int nshells = s_arg->shells.size()-1;
		assert(nshells >= 1);

        coordlist_t& x = shape.x;
        std::vector<double>& weight = shape.f;

        sd.resize(nsectors*nshells, component_t(0.0, 0.0));

		std::vector<double> rsq = s_arg->shells;
		for (unsigned int i=0; i<rsq.size(); i++) {
			rsq[i] *= rsq[i];
		}

		int ntotal = 0;
		for (int s=0; s<nshells; s++) {
			int n = x.size();
			shpdesc_t sd_shell(nsectors, component_t(0.0, 0.0));
			double twopi = 2.0*M_PI;

			for (int i=0; i<n; i++) {
				double rsqi = x[i][0]*x[i][0] + x[i][1]*x[i][1];

				if (rsqi >= rsq[s] && rsqi < rsq[s+1]) {
					double theta = atan2(x[i][1], x[i][0]) + M_PI;
					if (theta == twopi) {
                        theta -= twopi;
                    }
                    sd_shell[(int)(theta/twopi*nsectors)] += component_t(weight[i]*1.0);
					ntotal++;
				}
			}
            //std::cerr << sd.size() << "\t" << sd_shell.size() << std::endl;
			//sd.insert(sd.end(), sd_shell.begin(), sd_shell.end());
            for (int j=0; j<nsectors; j++) {
                sd[s*nshells + j] = sd_shell[j];
            }
		}

		if (s_arg->normalize) {
			for (unsigned int i=0; i<sd.size(); i++) {
				sd[i] /= (double)ntotal;
			}
		}
	}

	/**
	@param x is an array of coordinates
	@param nbin in the number of bins per one revolution around the sphere
	@param sd is an array of real numbers representing the computed shape
		descriptor
	*/
	void shphist3d(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		shphist_info* s_arg;

        shphist_info default_info;
        if (arg == 0x0) {
            s_arg = &default_info;
        }
        else {
            s_arg = static_cast<shphist_info*>(arg);
        }

        coordlist_t& x = shape.x;
        std::vector<double>& weight = shape.f;

		int nsectors = s_arg->nsectors;
		int nshells = s_arg->shells.size()-1;
		assert(nshells >= 1);

		std::vector<double> r_shell = s_arg->shells;

		int n=x.size(), ntotal=0;
		sd.clear();
		int nsectorsq = nsectors*nsectors;

		for (int s=0; s<nshells; s++) {
			shpdesc_t sd_shell(nsectorsq, 0.0);

			for (int i=0; i<n; i++) {

				double r = sqrt(x[i][0]*x[i][0] +
								x[i][1]*x[i][1] +
								x[i][2]*x[i][2]);

				if (r >= r_shell[s] && r < r_shell[s+1]) {
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

					int phi_bin = (int)( phi/M_PI*nsectors*0.5);							//phi [0, 2pi]
					int theta_bin = (int)(theta/M_PI*nsectors);							//theta [0,pi]

					sd_shell[theta_bin*nsectors + phi_bin] += component_t(weight[i]*1.0);
					ntotal++;
				}
                //sd.insert(sd.begin(), sd_shell.begin(), sd_shell.end());
            }
            sd.insert(sd.end(), sd_shell.begin(), sd_shell.end());
		}

		if (s_arg->normalize) {
			for (unsigned int i=0; i<sd.size(); i++) {
				sd[i] /= (double)ntotal;
			}
		}
	}

	void shphistico(coordlist_t& x, shpdesc_t& sd, arg_t arg)
	{
		//shphistico_info* s_arg = static_cast<shphistico_info*>(arg);

		/*
		if (s_arg->gridprecision != s_arg->oldgridprecision) {
			int factor = s_arg->gridprecision, node_num, edge_num, triangle_num;
			sphere_imp_grid_icos_size (factor, &node_num, &edge_num, &triangle_num);
			std::vector<double> xyz(node_num*3);
			sphere_imp_gridpoints_icos1 (factor, node_num, &xyz[0]);

			s_arg->cells.clear();
			s_arg->cells.setBox(2.0);
			s_arg->cells.setPeriodicity(
			s_arg->xnode.resize(node_num);
			for (int i=0; i<node_num; i++) {
				coord_t xi(3);
				for (int j=0; j<3; j++) {
					xi[j] = xyz[i*3+j];
				}
				s_arg->xnode[i] = xi;
				s_arg->cells.insert(i, xi);
			}
		}

		CellList& cells = s_arg->cells;
		for (unsigned int i=0; i<x.size(); i++) {
			coord_t xi = x[i];

			Cell& icell = cells.getCell(cells.getCellContaining(xi));

			s_arg->cells.getNeighborsOf(i);

		}
		*/
	}

    void shphistalign2d(shpdesc_t& sd1, shpdesc_t& sd2, shphist_info& info)
    {
        assert(sd1.size() == sd2.size());

		int nsectors = info.nsectors;
		int nshells = info.shells.size()-1;
		assert(nshells >= 1);

        shpdesc_t sdorig = sd1;
        shpdesc_t sdshift = sd1;
        double max_match = -10000.0;
        int size = sd1.size();
        for (unsigned int shift=0; shift<nsectors; shift++) {
            for (int s=0; s<nshells; s++) {
                for (unsigned int i=0; i<nsectors; i++) {
                    sdshift[s*nsectors + (i+shift)%nsectors] = sd1[s*nsectors + i];
                }
            }
            double match = matchfundist(sdshift, sd2);
            if (match > max_match) {
                sd1 = sdshift;
                max_match = match;
            }
        }
    }
}
