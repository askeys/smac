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
 *  centroid_distance.cpp
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "shpdescrms.h"
#include "preprocess.h"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "matchfun.h"

#ifdef _USE_HUNGARIAN
#include "../3rdparty/hungarian/hungarian.h"
#endif


#include <assert.h>
#include <cstdlib>


namespace smac
{
	/**
    \brief Computes the RMS descriptor
	\ingroup shpdesc
    \param shape is shapedata
    \param sd is a placeholder for the resulting shape descriptor
    \param arg is a pointer to a shpdescrms_info struct (or NULL)
    \note this descriptor is non-invariant and therefore requires size
        normalization and registration for scale/rotation-independent matching
    Computes the the RMS descriptor for a set of shape data.  The descriptor
    simply consists of the points concatenated into a descriptor.  If the weight
    of a point is greater than 1, the point is added multiple times.  Currently,
    the descriptor only supports integer weights.  The descriptor requires an
    assigment step before matching
    */
	void rmsdesc(shapedata_t& shape, shpdesc_t& sd, arg_t arg)
	{
		sd.clear();
		coordlist_t x = shape.x;
        weight_t w = shape.f;

        if (x.size() == 0) {
			return;
		}

        int dim = x[0].size();
		for (unsigned int i=0; i<x.size(); i++) {
            for (double t=0; t<w[i]; t+=1.0) {
                for (int j=0; j<3; j++) {
                    if (j >= dim) {
                        sd.push_back(component_t(0.0));
                    }
                    else {
                        sd.push_back(component_t(x[i][j], 0.0));
                    }
                }
			}
		}
	}

    /**
    \brief Figures out the assignment between two RMS descriptors and reorders
        one of the descriptors accordingly
    \ingroup shpdesc
    \param sd1 is an RMS descriptor
    \param sd2 is an RMS descriptor
    \param arg is a pointer to a rmsfastassign_info struct
    \bug currently only works for sd1.size() == sd2.size()

    The RMS fast assign method tries to assign points in the RMS descriptor
    based on their minimum distance or max dot product.  In general, minimizing
    the distance is optimal for scale-normalized shapes while maximizing the dot
    product is optimal for shapes of different sizes.  The assignment is
    non-unique; therefore, more than one point on shape1 may be assigned to a
    given point on shape2.  Unlike the more complex rmsassign method, this
    method has no problems with degeneracy.
    */
	void rmsfastassign(shpdesc_t& sd1, shpdesc_t& sd2, arg_t arg)
	{
        rmsfastassign_info* info = static_cast<rmsfastassign_info*> (arg);
        rmsfastassign_info default_info;
        if (arg == 0x0) {
            info = &default_info;
        }

        shpdesc_t sdtmp;

        int n = sd1.size()/3;
		int m = sd2.size()/3;

        //we always want m>=n
        //this way, we can loop over all n and assign an m
        bool swap = false;
        if (m<n) {
            sdtmp = sd1;
            sd1 = sd2;
            sd2 = sdtmp;
            int tmp = n;
            n = m;
            m = tmp;
            swap = true;
        }

        assert(m>=n);

        sdtmp.resize(m*3);
        std::vector<bool> counted(m, false);

		for (int i=0; i<n; i++) {
			double max_match = -1e10;
			int max_index = 0;
			for (int j=0; j<m; j++) {
                if (counted[j]) {
                    continue;
                }
                if (info->method == RMS || info->method == RSQ) {
                    double diff = 0.0;
                    for (int k=0; k<3; k++) {
                        double dx = sd1[i*3 + k].real() - sd2[j*3 + k].real();
                        diff += dx*dx;
                    }
                    if(-1.0*diff > max_match) {
                        max_index = j;
                        max_match = -1.0*diff;
                    }
                }
                else if (info->method == DOT) {
                    double dot = 0.0, m1 = 0.0, m2 = 0.0, norm;
                    for (int k=0; k<3; k++) {
                        dot += sd1[i*3 + k].real()*sd2[j*3 + k].real();
                        m1 += sd1[i*3 + k].real()*sd1[i*3 + k].real();
                        m2 += sd2[j*3 + k].real()*sd2[j*3 + k].real();
                    }
                    norm = m1*m2;
                    double dotp = dot / norm;
                    if(dotp > max_match) {
                        max_index = j;
                        max_match = dotp;
                    }
                }
                else {
                    std::cerr << "***ERROR: rmsfastassign method ";
                    std::cerr << "not recognized\n";
                    throw std::runtime_error("ERROR using rmsfastassign");
                }
			}
            if (counted[max_index]) {
                std::cerr << "***WARNING: bad descriptor(s) passed to ";
                std::cerr << "rmsfastassign; skipping assignment\n";
                return;
            }
            counted[max_index] = true;

            for (int k=0; k<3; k++) {
                sdtmp[i*3+k] = sd2[max_index*3+k];
            }
 		}
        int nc = n;
        for (int i=0; i<m; i++) {
            if (!counted[i]) {
                for (int k=0; k<3; k++) {
                    sdtmp[nc*3+k] = sd2[i*3+k];
                }
                nc++;
            }
        }

        if (swap) {
            sd2 = sd1;
            sd1 = sdtmp;
        }
        else {
            sd2 = sdtmp;
        }
	}

    /**
    \brief Figures out the assignment between two RMS descriptors and reorders
        one of the descriptors accordingly
    \ingroup shpdesc
    \param sd1 is an RMS descriptor
    \param sd2 is an RMS descriptor
    \param arg is a pointer to a rmsassign_info struct
    \bug hangs for degenerate shapes (where fitness is equal between points)

    The RMS assign method tries to assign points in the RMS descriptor
    based on their minimum distance or max dot product.  In general, minimizing
    the distance is optimal for scale-normalized shapes while maximizing the dot
    product is optimal for shapes of different sizes.  The function uses the
    Hungarian method to determine the optimum assignment.
    */
	void rmsassign(shpdesc_t& sd1, shpdesc_t& sd2, arg_t arg)
	{
#ifdef _USE_HUNGARIAN
        rmsfastassign_info* info = static_cast<rmsfastassign_info*> (arg);

		int n = sd1.size()/3;
		int m = sd2.size()/3;

		int *cost = new int [n*m];

		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				double dot = 0.0, m1 = 0.0, m2 = 0.0, norm;
				for (int k=0; k<3; k++) {
					dot += sd1[i*3 + k].real()*sd2[j*3 + k].real();
					m1 += sd1[i*3 + k].real()*sd1[i*3 + k].real();
					m2 += sd2[j*3 + k].real()*sd2[j*3 + k].real();
				}
				norm = m1*m2;

				if (fabs(norm) > 1e-2) {
					cost[i*n+j] = 10000*dot/norm;
				}
				else {
					cost[i*n+j] = 0;
				}
			}
		}

		hungarian_t prob;
		hungarian_init(&prob, cost, n, m, HUNGARIAN_MAX);
		if (info->verbose) {
			hungarian_print_rating(&prob);
		}
		hungarian_solve(&prob);

		for (int i=0; i<n; i++) {
			for (int j=i+1; j<m; j++) {
				if (j==prob.a[i] && i!=j) {
					component_t temp[3];
					temp[0] = sd2[i*3+0];
					temp[1] = sd2[i*3+1];
					temp[2] = sd2[i*3+2];
					sd2[i*3+0] = sd2[j*3+0];
					sd2[i*3+1] = sd2[j*3+1];
					sd2[i*3+2] = sd2[j*3+2];
					sd2[j*3+0] = temp[0];
					sd2[j*3+1] = temp[1];
					sd2[j*3+2] = temp[2];
				}
			}
		}

		if (info->verbose) {
			hungarian_print_assignment(&prob);
		}
		hungarian_fini(&prob);
		delete cost;
#else
        std::cerr << "***ERROR: rmsassign called, but Hungarian package was ";
        std::cerr << "disabled.  Try rmsfastassign instead.";
        throw std::runtime_error("Error calling Hungarian");
#endif
	}


	/**
	*/
	double matchfundotrms(shpdesc_t& s1, shpdesc_t& s2, double min, double max)
    {
        rmsfastassign(s1, s2);
        return matchfundot(s1, s2, min, max);
    }

	/**
	*/
	double matchfundistrms(shpdesc_t& s1, shpdesc_t& s2, double min, double max)
	{
        rmsfastassign(s1, s2);
        return matchfundist(s1, s2, min, max);
    }

	/**
	*/
	double matchfundiffrms(shpdesc_t& s1, shpdesc_t& s2, double min, double max)
    {
        rmsfastassign(s1, s2);
        return matchfundiff(s1, s2, min, max);
    }

	/**
	*/
	double matchfunchisqrms(shpdesc_t& s1, shpdesc_t& s2, double min, double max)
    {
        rmsfastassign(s1, s2);
        return matchfunchisq(s1, s2, min, max);
    }
}
