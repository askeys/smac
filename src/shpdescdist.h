/*
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
 *  shape_distribution.h
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SHAPEDISTRIBUTION
#define SHAPEDISTRIBUTION

#include "typedefs.h"

#ifdef __linux__
#include <assert.h>
#endif


namespace smac
{
    enum {COMPLEX=0, REAL=1, IMAG=2};
    
    struct shpdistD2_info
    {
        shpdistD2_info() : rmin(0.0), rmax(5.0), nbin(25) {}
        double rmin;
        double rmax;
        int nbin;
        box_t box;
    };

    struct shpdistA3_info
    {
        shpdistA3_info() : rmin(0.0), rmax(5.0), nbin(25) {}
        double rmin;
        double rmax;
        int nbin;
        box_t box;
    };

	void shpdistD2(shapedata_t& s, shpdesc_t& sd, arg_t);
	void shpdistA3(shapedata_t& s, shpdesc_t& sd, arg_t);
	void shpdist(shpdesclist_t&, shpdesc_t&, arg_t);

	struct shpdist_info
	{
		shpdist_info() :
			nbin(10), min(0.0), max(1.0), comptype(COMPLEX) {}

		int nbin;
        int comptype;                                                           ///<type of components (REAL, IMAG, COMPLEX)
		double min;
		double max;
	};
}

#endif
