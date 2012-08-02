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
 *  centroid_distance.h
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef RMSDESCRIPTOR
#define RMSDESCRIPTOR

#include "typedefs.h"

namespace smac
{
    /**
    \brief holds inputs for the RMS descriptor
    \ingroup shpdesc
    */
    struct rmsdesc_info 
    {
    };

    enum {RMS=1, RSQ=2, DOT=3}; 

    /**
    \brief holds inputs for the assignement step for the RMS descriptor
    \ingroup shpdesc
    */
    struct rmsassign_info 
    {
        rmsassign_info() : verbose(false), method(RMS) {}
        bool verbose;
        int method;
    };

    typedef rmsassign_info rmsfastassign_info;
    
	void rmsdesc(shapedata_t& x, shpdesc_t& sd, arg_t arg);                     
    void rmsassign(shpdesc_t& sd1, shpdesc_t& sd2, arg_t arg);                  
	void rmsfastassign(shpdesc_t& sd1, shpdesc_t& sd2, arg_t arg=0x0);            

    //specialized matching functions for RMS descriptors:
    
	/**
	\brief Aligns two RMS descriptors and computes the dot product using 
        matchfundot
	\ingroup matchfun
	*/
	double matchfundotrms(shpdesc_t&, shpdesc_t&,
		double min=0.0, double max=1.0);
		
	/**
	\brief Aligns two RMS descriptors and computes the euclidean distance 
        metric usigin matchfundist
	\ingroup matchfun
	*/
	double matchfundistrms(shpdesc_t&, shpdesc_t&, 
		double min=0.0, double max=1.0);
	
	/**
	\brief Aligns two RMS descriptors and computes the difference metric 
        using matchfundiff
    \ingroup matchfun
	*/
	double matchfundiffrms(shpdesc_t&, shpdesc_t&, 
		double min=0.0, double max=1.0);
	
	/**
	\brief Aligns two RMS descriptors and computes the chi-squared metric
        using matchfunchisq
	\ingroup matchfun
	*/
	double matchfunchisqrms(shpdesc_t&, shpdesc_t&,	
		double min=0.0, double max=1.0);	
}

#endif
