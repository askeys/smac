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

#include "rand.h"
#include <cmath>


#include <cstdlib>

namespace smac
{
    /**
    \brief generates a random number from a gaussian distribution
    \ingrou util
    \param mean is the mean of the distribution
    \param variance is the variance of the distribution
    */
	double randgauss(double mean, double variance)
	{
		double first,v1,v2,rsq,fac;
		static int save = 0;
		static double second;

		if (!save) {
			int again = 1;
			while (again) {
				v1 = 2.0*drand48()-1.0;
				v2 = 2.0*drand48()-1.0;
				rsq = v1*v1 + v2*v2;
				if (rsq < 1.0 && rsq != 0.0) again = 0;
			}
			fac = sqrt(-2.0*log(rsq)/rsq);
			second = v1*fac;
			first = v2*fac;
			save = 1;
		}
		else {
			first = second;
			save = 0;
		}
		return mean + first*variance;
	}

    /**
    \brief generates a random quaternion uniform in SO(3)
    \ingrou util
    \param q is a placeholder for the random quaternion
    */
    void randquaternion(std::vector<double>& q)
    {
        q.resize(4);
        double u1 = drand48();
        double u2 = drand48();
        double u3 = drand48();
        double twopiu2 = 2.0*M_PI*u2;
        double twopiu3 = 2.0*M_PI*u3;
        q[0] = sqrt(1.0-u1)*sin(twopiu2);
        q[1] = sqrt(1.0-u1)*cos(twopiu2);
        q[2] = sqrt(u1)*sin(twopiu3);
        q[3] = sqrt(u1)*cos(twopiu3);
    }
}
