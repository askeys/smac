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
 *  rdf.cpp
 *  libsmac
 *
 *  Created by askeys on 2/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "diagnostics.h"
#include "typedefs.h"
#include "space.h"
#include "HistogramDynamic.hpp"

#include <algorithm>


namespace smac
{
	/**
	\param x is the coordinates
	\param box is the box
	\param args is a struct of type rdf_info containing additional arguments
	\param gofr is the container to hold the output
	*/
	void rdf(coordlist_t& x, box_t& box, rdf_info& args, function1d_t& gofr)
	{
		double binsize = args.binsize;
		double rmax = args.rmax;

		HistogramDynamic <double> radial_histogram(binsize);
		radial_histogram.insert(0,0);											//initialize the histogram to always start at / line up with zero

		coord_t L(box.period);
		double Lmin = *(std::min_element(L.begin(), L.end()));
		double L2sq = Lmin * 0.5;
		L2sq *= L2sq;
		if (rmax > 0) {
			L2sq = rmax*rmax;
		}

		int natom = x.size();
		for (int i=0; i<natom; i++) {
			for (int j=i+1; j<natom; j++) {
				double rsq = distancesq(x[i], x[j], box.period, box.periodic);
				if (rsq < L2sq) {
					double r = sqrt(rsq);
					radial_histogram.insert(r, 2.0);							//give weight of 2, since we are doing 1/2 N^2 distances
				}
			}
		}

		double np = (double)natom;
		double delta = binsize / 2.0;
		double vol = 1;
		for (unsigned int i=0; i<box.period.size(); i++) {
			vol *= box.period[i];
		}
		double density = np / vol;

		function1d_t normalization, probability;
		for (unsigned int i=0; i<radial_histogram.getNumberOfBins(); i++) {
			std::pair<double, double> bin = radial_histogram.getBin(i);
			double r = bin.first;
			double count = bin.second;

			double vol_shell =	4.0/3.0*M_PI*pow(r+delta, 3) -
								4.0/3.0*M_PI*pow(r-delta, 3);
			double np_ideal_gas = density * vol_shell;

			if (normalization.find(r) == normalization.end()) {
				normalization[r] = 0.0;
				probability[r] = 0.0;
			}

			normalization[r] += np_ideal_gas * np;
			probability[r] += count;
		}

		for(function1d_t::iterator i=normalization.begin();
			i!=normalization.end(); ++i) {
			gofr[i->first] = probability[i->first] / i->second;
		}
	}
}
