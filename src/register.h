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
 *  register.h
 *  smac
 *
 *  Created by askeys on 12/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _SMAC_REGISTER_
#define _SMAC_REGISTER_

#include "typedefs.h"
#include "matchfun.h"
#include "shpdescfourier.h"

namespace smac
{
    enum ICPTrialType {LINEAR, RANDOM};
    /**
    \brief optional arguments for th registericp function
    \ingroup reg
    */
    struct registericp_info 
    {
        registericp_info() : maxiter(20), maxrotations(216), tol(0.0), 
            icptol(0.7), trialtype(LINEAR), matchfun(matchfundiff) {}
        int maxiter;                                                            ///<max number of ICP iterations
        int maxrotations;                                                       ///<max number of initial rotations
        double tol;                                                             ///<ICP stops if match > 1.0-tol
        double icptol;                                                          ///<Don't even attempt ICP if match < 1-tol
        int trialtype;
        matchfun_t matchfun;
    };
    
    void registericp(coordlist_t& x, coordlist_t& xref, arg_t=0x0);

    /**
    \brief optional arguments for th registericp function
    \ingroup reg    
    */
    struct registerpca_info
    {
        registerpca_info() : mode('r') {}
        char mode;
    };
    
	void registerpca(coordlist_t& x, arg_t=0x0);
    
    /**
    \brief optional arguments for the registersymmetryaxis function
    \ingroup reg
    */
    struct registersymmetryaxis_info
    {
        registersymmetryaxis_info() : niter(30), nstartpoints(1000), 
            dqmax(0.02), beta(20), eterminate(-1e20), fourierargs(0x0) {}
        
        int niter;
        int nstartpoints;
        double dqmax;
        double beta;
        double eterminate;
        fourierdesc2d_info* fourierargs;
    };

    void registersymmetryaxis(
        coordlist_t& x_orig, 
        std::vector<double>& q_best, 
        registersymmetryaxis_info& info);
    
    double registerzerocomplexcomonent(coordlist_t& x_orig, fourier_info& info);


}

#endif
