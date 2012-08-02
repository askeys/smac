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
 *  matching_metric.h
 *  smac
 *
 *  Created by askeys on 12/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATCHINGMETRIC
#define MATCHINGMETRIC

#include "typedefs.h"

namespace smac 
{
	/**
	\brief matching function based on the dot product of one shape descriptor
		with another
	\ingroup matchfun
	*/
	double matchfundot(shpdesc_t&, shpdesc_t&,
		double min=0.0, double max=1.0);
		
	/**
	\brief matching function based on the euclidean distance between shape
		descriptors
	\ingroup matchfun
	*/
	double matchfundist(shpdesc_t&, shpdesc_t&, 
		double min=0.0, double max=1.0);
	
	/**
	\brief matching function based on the absolute difference between shape
		descriptors
	\ingroup matchfun
	*/
	double matchfundiff(shpdesc_t&, shpdesc_t&, 
		double min=0.0, double max=1.0);
	
	/**
	\brief matching function based on the chi-squared metric for two shape
		descriptors
	\ingroup matchfun
	*/
	double matchfunchisq(shpdesc_t&, shpdesc_t&,	
		double min=0.0, double max=1.0);
	
	struct matchvalue_t
	{
		int index;
		double match;
	};
	
	typedef std::vector<matchvalue_t> matchlist_t;
	void bestmatch(shpdesc_t&, shpdesclist_t&, matchfun_t, matchlist_t&, 
        assignfun_t=0x0, arg_t=0x0);

}

#endif
