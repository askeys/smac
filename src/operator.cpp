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
 *  operations.cpp
 *  smac
 *
 *  Created by askeys on 12/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "operator.h"
#include <cmath>

namespace smac {
	

    void divide(std::vector<double>& v, double d)
    {
        for (unsigned int i=0; i<v.size(); i++) {
            v[i] /= d;
        }
    }
    void multiply(std::vector<double>& v, double d)
    {
        for (unsigned int i=0; i<v.size(); i++) {
            v[i] *= d;
        }
    }

    void plusequals(std::vector<double>& v1, std::vector<double>& v2)
    {
        for (unsigned int i=0; i<v1.size(); i++) {
            v1[i] += v2[i];
        }
    }

    void minusequals(std::vector<double>& v1, std::vector<double>& v2)
    {
        for (unsigned int i=0; i<v1.size(); i++) {
            v1[i] -= v2[i];
        }
    }
    

    void normalize(std::vector<double>& v)
    {
        double sum = 0.0;
        for (unsigned int i=0; i<v.size(); i++) {
            sum += v[i]*v[i];
        }
        sum = sqrt(sum);
        for (unsigned int i=0; i<v.size(); i++) {
            v[i] /= sum;
        }
    }

	/**
    \brief Squared magnitude of a shape descriptor vector
    \ingroup util
	*/
	double magnitudesq(shpdesc_t& sd)
	{
		double sum_sq = 0.0;
		for (unsigned int i=0; i<sd.size(); i++) {
			sum_sq += sd[i].real()*sd[i].real() 
							+ sd[i].imag()*sd[i].imag();
		}
		return sum_sq;
	}

	/**
    \brief Magnitude of a shape descriptor vector
    \ingroup util
	*/
	double magnitude(shpdesc_t& sd)
	{
		return sqrt(magnitudesq(sd));
	}

	/**
    \brief Add two shape descriptor vectors
    \ingroup util
	*/
	void add(shpdesc_t& sd1, shpdesc_t& sd2) 
	{
		sd1.resize(std::max(sd1.size(), sd2.size()));
		
		for (unsigned int i=0; i<sd2.size(); i++) {
			sd1[i] += sd2[i];
		}
	}

	/**
    \brief Subtract two shape descriptor vectors
    \ingroup util
	*/
	void subtract(shpdesc_t& sd1, shpdesc_t& sd2) 
	{
		sd1.resize(std::max(sd1.size(), sd2.size()));
		
		for (unsigned int i=0; i<sd2.size(); i++) {
			sd1[i] -= sd2[i];
		}
	}
	
	/**
    \brief Divide a shape descriptor vector by a scalar
    \ingroup util
	*/
	void divide(shpdesc_t& sd, double d)
	{
		for (unsigned int i=0; i<sd.size(); i++) {
			sd[i] /= d;
		}		
	}

	/**
    \brief Multiply a shape descriptor vector by a scalar
    \ingroup util
	*/
	void multiply(shpdesc_t& sd, double d)
	{
		for (unsigned int i=0; i<sd.size(); i++) {
			sd[i] *= d;
		}		
	}
	
	/**
    \brief Normalize a shape descriptor vector to unit length
    \ingroup util
	*/
	void normalize(shpdesc_t& sd)
	{
	
		double sum = 0.0;
		for (unsigned int i=0; i<sd.size(); i++) {
			sum += sd[i].real()*sd[i].real() + sd[i].imag()*sd[i].imag();
		}
		
		for (unsigned int i=0; i<sd.size(); i++) {
			sd[i] /= sqrt(sum);
		}
	}

	/**
    \brief Concatenate shape descriptor vectors
    \ingroup util
	*/
	void aggregate(shpdesc_t& sd1, shpdesc_t& sd2) 
	{
		sd1.insert(sd1.end(), sd2.begin(), sd2.end());
	}
}
