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
 *  matching_metric.cpp
 *  smac
 *
 *  Created by askeys on 12/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "matchfun.h"
#include "operator.h"
#include <iostream>


#include <assert.h>
#include <cstdlib>

namespace smac
{
	/**
	\brief Matches descriptors based on the euclidean distance metric
	\ingroup matchfun
	Normalized Sqrt of sum of square differences, similar to chi^2, R norm
	*/
	double matchfundist(
		shpdesc_t& sd1, shpdesc_t& sd2,
		double min, double max)
	{
		assert(sd1.size() == sd2.size());

		double max_difference = magnitude(sd1) + magnitude(sd2);
		if (max_difference == 0.0) {
			return max;
		}

		double residual = 0.0;
		std::complex<double> diff_i;

		for (unsigned int i=0; i<sd1.size(); i++) {
			diff_i = sd1[i] - sd2[i];
			residual += diff_i.real() * diff_i.real()
					  + diff_i.imag() * diff_i.imag();
		}
		double r = sqrt(residual);

		double rnorm = r / max_difference;
		double range = max - min;
		double rmm = min + range*(1.0-rnorm);
		return rmm;
	}

	/**
	\brief Matches descriptors based on diff metric
	\ingroup matchfun
	Normalized sum of magnitude of differences, D_{-1,1}
	*/
	double matchfundiff(
		shpdesc_t& sd1, shpdesc_t& sd2,
		double min, double max)
	{
		assert(sd1.size() == sd2.size());

		double den = 0.0;
		for (unsigned int i=0; i<sd1.size(); i++) {
			den += std::abs(sd1[i]) + std::abs(sd2[i]);
		}
		if (den == 0.0) {
            std::cerr << "***WARNING: matchfundiff: zero descriptor\n";
			return max;
		}

		double residual = 0.0;
		std::complex<double> diff_i;

		for (unsigned int i=0; i<sd1.size(); i++) {
			diff_i = sd1[i] - sd2[i];
			residual += std::abs(diff_i);
		}

		double range = max - min;
		double dnorm =  min + range * (1.0 - residual/den);
		return dnorm;
	}

	/**
	\brief Matches descriptors based on the scalar product metric
	\ingroup matchfun
	Normalized dot product, P_{-1,1}
	*/
	double matchfundot(
		shpdesc_t& sd1, shpdesc_t& sd2,
		double min, double max)
	{
		assert(sd1.size() == sd2.size());

		double max_p = magnitude(sd1) * magnitude(sd2);

		double dot_product = 0.0;
		for (unsigned int i=0; i<sd1.size(); i++) {
			dot_product += sd1[i].real()*sd2[i].real()
							+ sd1[i].imag()*sd2[i].imag();
		}

		double range = max - min;
		double dot_norm = range*(dot_product/max_p + 1.0)*0.5 + min;
		return dot_norm;
	}

	/**
	\brief Matches descriptors based on the chisq metric
	\ingroup matchfun
	0.5 x the normalized sum of the squares, chi^2.  Also sometimes refered to
	as the cost function
	*/
	double matchfunchisq(
		shpdesc_t& sd1, shpdesc_t& sd2,
		double min, double max)
	{
		assert(sd1.size() == sd2.size());

		double max_difference = 0.0;
		for (unsigned int i=0; i<sd1.size(); i++) {
            double tmp = std::abs(sd1[i]) + std::abs(sd2[i]);
            max_difference += tmp*tmp;
        }

        if (max_difference == 0.0) {
            std::cerr << "***WARNING: matchfunchisq: zero descriptor\n";
			return max;
		}

		double residual = 0.0;
		std::complex<double> diff_i;

		for (unsigned int i=0; i<sd1.size(); i++) {
			diff_i = sd1[i] - sd2[i];
			residual += diff_i.real() * diff_i.real()
					  + diff_i.imag() * diff_i.imag();
		}

		double rnorm = residual / max_difference;
		double range = max - min;
		double rmm = min + range*(1.0-rnorm);
		return rmm;
	}

	/**
	\brief Finds the best match for a list of reference structures
	\ingroup matchfun
    \param s is the shape descriptor of the unknown
    \param sref is a list of reference descritpors
    \param matchfun is the matching function to use
	\param match is an ordered list to hold the index / matching value pairs
    \param assignfun is the assignment function to use, if applicable
    \param assignarg is the a pointer to the argument struct associated with the
        assign function
	*/
	void bestmatch(shpdesc_t& s, shpdesclist_t& sref, matchfun_t matchfun,
        matchlist_t& match, assignfun_t assignfun, arg_t assignarg)
	{
		std::map<double, int> m;
		for (unsigned int i=0; i<sref.size(); i++) {
            if (assignfun) {
                assignfun(s, sref[i], assignarg);
            }
            m[matchfun(s, sref[i], 0.0, 1.0)] = i;
		}

		match.clear();
		for (std::map<double, int>::reverse_iterator i=m.rbegin(); i!=m.rend(); ++i) {
			matchvalue_t mi;
			mi.index = i->second;
			mi.match = i->first;
			match.push_back(mi);
		}
	}
}
