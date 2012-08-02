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
 *  MoveVelocityVerlet.h
 *  glotzilla
 *
 *  Created by Aaron Keys on 7/13/08.
 *  Copyright 2008 __MyCompanynatomame__. All rights reserved.
 *
 */

#ifndef SMCLUSTER
#define SMCLUSTER

#include "typedefs.h"
#include <vector>
#include <map>

namespace smac 
{		
	
	void clstrshortrange(coordlist_t& x, box_t&, double interaction_range, 
		clstrrule_t, arg_t args, clstrlist_t& cluster);
		
	void clstr(int n, clstrrule_t, arg_t args, clstrlist_t&);
		
	void clstrhierarchical(int n, clstrrule_t, arg_t args, clstrlist_t&);

	void clstrshortrangehierarchical(
		coordlist_t& x, box_t& box, 
		double interaction_range, clstrrule_t, arg_t, clstrlist_t&);


	//clustering rules
	bool clstrrulerangetype(int i, int j, arg_t);
	bool clstrrulerange(int i, int j, arg_t);
	bool clstrruleshell(int i, int j, arg_t);

	/**
	\brief argument struct for the clstrrulerange cluster function
	*/
	struct clstrrulerangetype_info {
		coordlist_t x;
		box_t box;
		unsigned int dim;
        std::vector<unsigned int> type;
		std::vector< std::vector<double> > rcutsq;

		clstrrulerangetype_info() : 
			dim(3)
			{}			
	};

	/**
	\brief argument struct for the clstrrulerange cluster function
	*/
	struct clstrrulerange_info {
		coordlist_t x;
		box_t box;
		int dim;
		double rcutsq;

		clstrrulerange_info() : 
			dim(3),
			rcutsq(0.0)
			{}			
	};

	/**
	\brief argument struct for the clstrruleshell cluster function
	*/
	struct clstrruleshell_info {
		coordlist_t x;
		box_t box;
		int dim;
		double rcutsqinner;
		double rcutsqouter;

		clstrruleshell_info() : 
			dim(3),
			rcutsqinner(0.0),
			rcutsqouter(0.0)
			{}			
	};
}

#endif
