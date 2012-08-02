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
 *  test_sd.cpp
 *  glotzilla
 *
 *  Created by askeys on 9/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 *  To compile: g++ test_sd.cpp -Iglotzilla_path glotz_lib_path/libglotzilla.a 
 */

#include <iostream>
#include <sstream>
#include <smac/smac.h>

using namespace std;
using namespace smac;

int main(int argc, char** argv)
{		
    string file;
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/helix/";
    cerr << datapath << endl;
    //we can watch the movie in the data directory to see that the last
    //datapoint is well formed (alternatively we could mathematically construct
    //an ideal helix) 
    xyzfile xyz;
    file = datapath + "backbone_51000000.xyz";
    loadxyz(file.c_str(), xyz);
    coordlist_t x = xyz.x;

	zernike_info args;
    args.rmax = 0.9;
	args.moment.clear();
	for (int i=4; i<=9; i++) {
		args.moment.push_back(i);
	}
	
	shpdesc_t descriptor, descriptor_ref;											
    shapedata_t shape_ref(x);
	zerndesc3d(shape_ref, descriptor_ref, &args);

	for (int i=1000000; i<51000000; i+=1000000) {
		ostringstream filename;
		filename << datapath << "backbone_" << i << ".xyz";
		loadxyz(filename.str().c_str(), xyz);
        x = xyz.x;
        shapedata_t shape_unknown(x);
		zerndesc3d(shape_unknown, descriptor, &args);

		//cerr << descriptor_ref.size() << endl;
		//cerr << descriptor.size() << endl;
        
		cout << i <<"\t"<< matchfundist(descriptor, descriptor_ref) << endl;
	}
	
	return 0;
}
