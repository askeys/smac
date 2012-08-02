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
 *  loadxyz_test.h
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
    cout << "XYZ reader Test\n";
    cout << "---------------\n";

    int errors = 0;
    
    ofstream testfile("test.xyz");
    testfile << "100\n";
    testfile << "10 10 10\n";
    for (unsigned int i=0; i<100; i++) {
        testfile << "H\t" << i << "\t" << i << "\t" << i << "\n";
    }
    testfile.close();
    
    xyzfile xyzdata;
    loadxyz("test.xyz", xyzdata);
    
    if (xyzdata.x.size() != 100) {
        cerr << "***ERROR: incorrect particle count for xyz file\n"; 
        errors++;
    }
    if (xyzdata.commentstr != "10 10 10") {
        cerr << "***ERROR: incorrect comment string for xyz file\n"; 
        errors++;
    }
    for (int i=0; i<100; i++) {
        if (xyzdata.type[i] != "H") {
            cerr << "***ERROR: incorrect type for xyz file ";
            cerr << "particle " << i << "\n"; 
            errors++;
        }
    }
    for (int i=0; i<100; i++) {
        for (int k=0; k<3; k++) {
            if ((int)xyzdata.x[i][k] != i) {
                cerr << "***ERROR: incorrect particle position for xyz file ";
                cerr << "particle " << i << "\n"; 
                errors++;
            }
        }
    }
    
    //system("rm test.xyz");
    
    if (errors) {
        cout << errors << " error(s) detected.\n";        
    }
    else {
        cout << "No errors detected.\n";
    }
    return 0;
}