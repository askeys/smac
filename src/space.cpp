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
 *  f_space.cpp
 *  glotzilla
 *
 *  Created by askeys on 7/16/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "space.h"
#include <cmath>

namespace smac 
{
	double distancesq(
		coord_t& x1, coord_t& x2, box_t& box)
	{
		int dim = x1.size();
		coord_t dx(dim);
		
		double sum = 0.0;
		for (int k=0; k<dim; k++) {
			dx[k] = x1[k] - x2[k];
			pbc(dx[k], box.period[k], box.periodic[k]);
			sum += dx[k]*dx[k];
		}
		return sum;
	}


	double distancesq(
		coord_t& x1, coord_t& x2, coord_t& period, bool_array_t& periodic)
	{
		int dim = x1.size();
		coord_t dx(dim);
		
		double sum = 0.0;
		for (int k=0; k<dim; k++) {
			dx[k] = x1[k] - x2[k];
			pbc(dx[k], period[k], periodic[k]);
			sum += dx[k]*dx[k];
		}
		return sum;
	}

	coord_t distancevec(coord_t& x1, coord_t& x2, 
		coord_t& period, bool_array_t& periodic)
	{
		int dim = x1.size();
		coord_t dx(dim);
		
		for (int k=0; k<dim; k++) {
			dx[k] = x1[k] - x2[k];
			pbc(dx[k], period[k], periodic[k]);
		}
		return dx;
	}

	double distance(coord_t& x1, coord_t& x2, 
		coord_t& period, bool_array_t& periodic)
	{
		return sqrt(distancesq(x1, x2, period, periodic));
	}

	double distancesq(coord_t& x1, coord_t& x2)
	{
		int dim = x1.size();
		coord_t dx(dim);
		
		double sum = 0.0;
		for (int k=0; k<dim; k++) {
			dx[k] = x1[k] - x2[k];
			sum += dx[k]*dx[k];
		}
		return sum;
	} 

	double distance(coord_t& x1, coord_t& x2)
	{
		return sqrt(distancesq(x1, x2));
	}

	void pbc(coord_t& x, coord_t& period, bool_array_t& periodic)
	{
		int dim = x.size();
		for (int k=0; k<dim; k++) {
			pbc(x[k], period[k], periodic[k]);
		}
	}

	inline double anint(double x)
	{
		if (x > 0.5) {
			return 1.0;
		}
		else if (x < -0.5) {
			return -1.0;
		}
		else {
			 return 0.0;
		}
	}
	
	void pbc(double& x, double period, bool periodic)
	{
		if (periodic) {
			x -= period*anint(x/period);
		}
	}

    bool inbox(coord_t& x, box_t box)
    {
        if (x[0] < box.boxlo[0]) {
            return false;
        }
        if (x[0] >= box.boxhi[0]) {
            return false;
        }
        if (x[1] < box.boxlo[1]) {
            return false;
        }
        if (x[1] >= box.boxhi[1]) {
            return false;
        }
        if (x[2] < box.boxlo[2]) {
            return false;
        }
        if (x[2] >= box.boxhi[2]) {
            return false;
        }
        return true;
    }
    
    bool inbox(coordlist_t& x, box_t box)
    {
        for (unsigned int i=0; i<x.size(); i++) {
            if (!inbox(x[i], box)) {
                return false;
            }
        }
        return true;
    }

}
