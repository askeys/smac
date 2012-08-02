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
 *  registericp_test.h
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

int main(int argc, char** argv)
{	
    cout << "ICP Registration Test\n";
    cout << "---------------\n";

    int errors = 0;
    
    coordlist_t x, xrotated;
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/";
    string file = datapath + "fcc_13.txt";
    load(file.c_str(), x);
    
    xrotated = x;
    rotate(xrotated, M_PI*drand48()*.5, M_PI*drand48()*-0.3, M_PI*drand48()*.3);
    
    shpdesc_t descriptor, descriptorrotated;
    shapedata_t shape(x), shaperotated(xrotated);
    
    rmsdesc_info sdinfo;
    rmsdesc(shape, descriptor, &sdinfo);
    rmsdesc(shaperotated, descriptorrotated, &sdinfo);
        
    double tol = 0.025;
    if (matchfundiff(descriptor, descriptorrotated, 0.0, 1.0) > 1.0-tol) {
        cerr << "***WARNING: descriptors optimally correspond before ";
        cerr << "ICP registration.\n";
        cerr << "(match value = " << matchfundiff(descriptor, descriptorrotated, 0.0, 1.0) << ")\n";
        errors++;    
    }
     
    registericp_info info;
    info.maxiter = 20;
    info.maxrotations = 200;
    info.trialtype = RANDOM;
    info.tol = tol;
    info.matchfun = matchfundiff;
    registericp(xrotated, x, &info);
    shaperotated.x = xrotated;
    rmsdesc(shaperotated, descriptorrotated, &sdinfo);
    rmsassign_info assigninfo;
    rmsfastassign(descriptor, descriptorrotated, &assigninfo);    
    
    if (matchfundiff(descriptor, descriptorrotated, 0.0, 1.0) < 1.0-tol) {
        cerr << "***ERROR: ICP registration failed ";
        cerr << "(match value = " << matchfundiff(descriptor, descriptorrotated, 0.0, 1.0) << ")\n";
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