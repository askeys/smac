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
 *  analysis.h
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef STANDARDANALYSIS
#define STANDARDANALYSIS

#include "typedefs.h"

namespace smac
{
    enum {CENTROID=1, CENTER_PARTICLE=2, SKIP_TRANSLATE=3};
    enum {MATCHSTD=0, MATCHNEGATE=1, MATCHABS=2};
    /**
    \brief holds optional arguments for the localorder function
    */
    struct localorder_info
    {
        localorder_info() : origin(CENTROID) {}
        int origin;
    };
    
	void localorder(coordlist_t& x, box_t& box, double range, 
		shpdescfun_t descfun, arg_t descargs, clstrlist_t&, shpdesclist_t&,
        localorder_info& info);
		
	void localorder_central_atom(coordlist_t& x, box_t& box, double range, 
		shpdescfun_t descfun, arg_t descargs, clstrlist_t& cluster,
		shpdesclist_t& descriptors);
	
	void bod_central_atom(coordlist_t& x,  box_t& box, double range, 
		shpdescfun_t descfun, arg_t descargs, clstrlist_t& cluster,
		shpdesc_t& descriptor, coordlist_t& x_combined);	
				
	struct crystalanalysis_info {
		crystalanalysis_info() :
		savedescriptors(false),
		saveclusters(false),
        matchmode(MATCHSTD) {}

		double rangedescriptor;													///<the range over which to calculate the descriptor
		double rangecompare;													///<the range over which to compare descriptors
		double rangecluster;													///<the range over which particles are considered to belong the same grain
		shpdescfun_t descriptor;												///<the shape descriptor function
		arg_t descriptorarg;													///<arguments for the shape descriptor function
		matchfun_t matchfun;													///<the matching function 
		double crystalbondcut;													///<cutoff for determining a a crystal-like bond
		double crystalparticlecut;												///<cutoff for the number of crystal-like bonds that defines a crystal particle 
		bool savedescriptors;													///<whether to save the descriptors used for the analysis
		bool saveclusters;														///<whether to save the clusters used to compute the descriptors
		shpdesclist_t saveddescriptors;											///<will contain the descriptors for each particle, if savedescriptors = true 
		clstrlist_t savedclusters;												///<will contain the clusters for each particle, if savedeclusters = true
        int matchmode;                                                          ///<whether to invert the match (for crystals whose neighbors mirror each other)
	};

	void testcrystalcriterion(coordlist_t& x, box_t& box, 
		double range1, double range2, shpdescfun_t descfun, 
		arg_t descargs, matchfun_t matchfun,
		double bondcut, function1d_t& dothistf, function1d_t& bondhistf,
        double histstep=0.05);

	void crystalorliquid(coordlist_t& x, box_t& box, 
		double range1, double range2,
		shpdescfun_t descfun, arg_t descargs, matchfun_t matchfun, 
		double bondcut, double fbondcut, 
		std::vector<bool>&);

	void crystalgrains(coordlist_t& x, box_t& box, 
		crystalanalysis_info&, clstrlist_t& grains);

}

#endif
