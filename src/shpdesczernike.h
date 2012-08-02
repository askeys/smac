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
 *  zerndesc.h
 *  libsmac
 *
 *  Created by askeys on 2/5/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ZERNDESC
#define ZERNDESC

#include "typedefs.h"

namespace smac
{
	void zerndesc2d(shapedata_t& x, shpdesc_t& sd, arg_t arg);
	void zerndesc3d(shapedata_t& x, shpdesc_t& sd, arg_t arg);
	void zernmom2d(shapedata_t& shape, shpdesc_t& sd, arg_t arg);
	void zernmom3d(shapedata_t& shape, shpdesc_t& sd, arg_t arg);

    enum {RADIALNORMCENTROID=1, 
        RADIALEXCLUDECENTROID, RADIALNORM, RADIALEXCLUDE, NONE}; 

	/**
	\brief A struct to hold information for computing Zernike moments
	*/
	struct zernike_info
	{
		zernike_info() : 
		invariant(true),
        normmethod(RADIALNORMCENTROID),
        rmax(0.99),
		_n(0),
		_l(0),
		_m(0)
		{
			moment.resize(11);
			for(int i=0; i<=10; i++) {
				moment[i] = i;
			}
		}
		
		std::vector<int> moment;                                                ///<List of Zernike moments to compute
		bool invariant;                                                         ///<Whether or not to compute invariant moments
        int normmethod;                                                         ///<How to normalize the shape
        double rmax;                                                            ///<Max radius of normalized shape
		int _n, _l, _m;															///<Pseudo-private variables for use with zermom functions
	};
    
    typedef zernike_info zerndesc2d_info;
    typedef zernike_info zerndesc3d_info;
}

#endif
