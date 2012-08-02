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
 *  fourierdescriptor_test.h
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include "smac.h"

using namespace std;
using namespace smac;

bool similar(double d1, double d2, double tol) 
{
    if (fabs(d1-d2) < tol) {
        return true;
    }
    else {
        return false;
    }
}

int main(int argc, char** argv)
{	
    cout << "Fourier Descriptor Test\n";
    cout << "---------------\n";

    int errors = 0;

    coordlist_t x_unknown, x_fcc, x_ico, x_hcp;
    string file;
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/basic_clusters/";    
    file = datapath + "fcc_13.txt";
    load(file.c_str(), x_fcc);
    file = datapath + "hcp_13.txt";
    load(file.c_str(), x_hcp);
    file = datapath + "ico_13.txt";
    load(file.c_str(), x_ico);
	
    shapedata_t shape_fcc(x_fcc), shape_hcp(x_hcp), shape_ico(x_ico);
    
	fourier_info args;
	args.frequency = vector<int>(1, 6);
	args.invariant = INVARIANT_Q;

    shpdesc_t q6_fcc, q6_hcp, q6_ico;
    
    //test values:
    fourierdesc3d(shape_fcc, q6_fcc, &args);
    fourierdesc3d(shape_hcp, q6_hcp, &args);
    fourierdesc3d(shape_ico, q6_ico, &args);
    
    if (!similar(q6_hcp[0].real(), 0.484, 0.01)) {
        cerr << "***ERROR: bad value for q6 hcp: " << q6_hcp[0] << "\n";
        errors++;    
    }
    if (!similar(q6_fcc[0].real(), 0.574524, 0.01)) {
        cerr << "***ERROR: bad value for q6 fcc: " << q6_fcc[0] << "\n";
        errors++;    
    }
    if (!similar(q6_ico[0].real(), 0.66, 0.01)) {
        cerr << "***ERROR: bad value for q6 ico: " << q6_ico[0] << "\n";
        errors++;    
    }
    
	//test rotational invariance:
	for (int i=0; i<5; i++) {
		rotate(x_ico, drand48(), drand48(), drand48());
        shapedata_t shape_ico_rot(x_ico);
        shpdesc_t descriptor_rot;
		fouriercoeff3d(shape_ico_rot, descriptor_rot, &args);
        if (!similar(q6_ico[0].real(), 0.66, 0.01)) {
            cerr << "***ERROR: invariance failed; q6 ico: " << q6_ico[0] << "\n";
            errors++;    
        }
	}
    
    //test shells:
    coordlist_t x_outer, x_inner;
    x_outer = x_ico;
    x_inner = x_ico;
	rescale(x_outer, 2.0);	
    vector<double> shells(3);
	shells[0] = 0.1; shells[1] = 1.2; shells[2] = 20;
	args.shells = shells;
	
    coordlist_t x_shells;
	x_shells.insert(x_shells.end(), x_inner.begin(), x_inner.end());
	x_shells.insert(x_shells.end(), x_outer.begin(), x_outer.end());
	
	shapedata_t shape_shells(x_shells);
    shpdesc_t q6_shell;
    fouriercoeff3d(shape_shells, q6_shell, &args);
    
    if (!similar(q6_shell[0].real(), 0.66, 0.01)) {
        cerr << "***ERROR: bad value for inner shell: " << q6_shell[0] << "\n";
        errors++;    
    }

    if (!similar(q6_shell[1].real(), 0.66, 0.01)) {
        cerr << "***ERROR: bad value for inner shell: " << q6_shell[1] << "\n";
        errors++;    
    }

    if (errors) {
        cout << errors << " error(s) detected.\n";        
    }
    else {
        cout << "No errors detected.\n";
    }
    return 0;
}
