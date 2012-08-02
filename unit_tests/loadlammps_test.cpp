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
 *  loadlammps_test.h
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
    cout << "LAMMPS dump reader Test\n";
    cout << "---------------\n";

    int errors = 0;
    
    ofstream testfile("test.lammpstrj");
    testfile << "ITEM: TIMESTEP\n";
    testfile << "1234\n";
    testfile << "ITEM: NUMBER OF ATOMS\n";
    testfile << "100\n";
    testfile << "ITEM: BOX BOUNDS\n";
    testfile << "0 100\n";
    testfile << "0 100\n";
    testfile << "0 100\n";
    testfile << "ITEM: ATOMS id type xs ys zs\n";

    for (unsigned int i=0; i<100; i++) {
        testfile << i+1 << "\t" << ((i<50)?1:2) << "\t" << i/100.0 << "\t" << i/100.0 << "\t" << i/100.0 << "\n";
    }
    testfile.close();
    
    lammpstrjfile lammpstrjdata;
    loadlammpstrj("test.lammpstrj", lammpstrjdata);

    if (lammpstrjdata.n != 100) {
        cerr << "***ERROR: incorrect n particles for lammpstrj file\n"; 
        errors++;
    }
    if (lammpstrjdata.x.size() != 100) {
        cerr << "***ERROR: incorrect particle count for lammpstrj file\n"; 
        errors++;
    }
    if (lammpstrjdata.timestep != 1234) {
        cerr << "***ERROR: incorrect comment string for lammpstrj file\n"; 
        errors++;
    }
    bool any_error = false;
    double tol = 0.01;
    for (int i=0; i<100; i++) {
        for (int k=0; k<3; k++) {
            if (i <= 50) {
                if (fabs(lammpstrjdata.x[i][k] - i) > tol) {
                    cerr << "***ERROR: incorrect particle position for lammpstrj file ";
                    cerr << "particle " << i << "\n"; 
                    any_error = true;
                }
            }
            else {
                if (fabs(lammpstrjdata.x[i][k] - (i-100)) > tol) {
                    cerr << "***ERROR: incorrect particle position for lammpstrj file ";
                    cerr << "particle " << i << "\n"; 
                    any_error = true;
                }
            }
        }
    }
    if (any_error) {
        errors++;
    }

    any_error = false;
    for (int i=0; i<100; i++) {
        if (lammpstrjdata.type[i] != ((i<50)?1:2)) {
            cerr << "***ERROR: incorrect type for lammpstrj file ";
            cerr << "particle " << i << "\n"; 
            any_error = true;
        }
    }
    if (any_error) {
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