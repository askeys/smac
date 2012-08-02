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
 *  visualization.cpp
 *  smac
 *
 *  Created by askeys on 12/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "neighbormap.h"

#include "space.h"
#include "io.h"
#include <cassert>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>

#include "cluster.h"
#include "preprocess.h"

namespace smac
{	
    /**
    \param shape is the input data
    \param icut is the inner radial cutoff
    \param ocot is the outer radial cutoff
    \param box is the boundary
    \param map is the output map
    */
    void neighbormap(
        shapedata_t& shape, double icut, double ocut, 
        box_t& box, shapedata_t& map)
	{
        coordlist_t& x = shape.x;
        map.x.clear();
        map.f.clear();
		assert(x.size() > 0);		
		clstrlist_t cluster;
		clstrruleshell_info args;
		args.x = x;
		args.box = box;
		args.dim = 3;
		args.rcutsqinner = icut*icut;
		args.rcutsqouter = ocut*ocut;
		
		cluster.clear();
		clstrshortrange(x, box, ocut, clstrruleshell, &args, cluster);
		
		for (unsigned int i=0; i<cluster.size(); i++) {
			if (cluster[i].size() <= 0) {
				continue;
			}

			coordlist_t xcluster(cluster[i].size()-1);
			
			for (unsigned int j=1; j<cluster[i].size(); j++) {
				xcluster[j-1] = x[cluster[i][j]];
			}			
            coord_t xc = x[cluster[i][0]], origin(3, 0.0);
            xc[0] = -xc[0];
            xc[1] = -xc[1];
            xc[2] = -xc[2];
            translate(xcluster, xc);
            unmap(xcluster, origin, box);
            //normusphere(xcluster);
                        
            for (unsigned ii=0; ii<xcluster.size(); ii++) {
                map.x.push_back(xcluster[ii]);
            }
			for (unsigned int j=0; j<cluster[i].size(); j++) {
                int ii = cluster[i][j];
                map.f.push_back(shape.f[ii]);
            }
		}
    }
}
