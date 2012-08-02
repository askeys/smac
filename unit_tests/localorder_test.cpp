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
 *  localorder_test.h
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
    cout << "Local Order Algorithm Test\n";
    cout << "---------------\n";

    int errors = 0;

    //set up particles on a cubic lattice
    coordlist_t x;
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            for (int k=0; k<10; k++) {
                coord_t xi(3);
                xi[0] = i;
                xi[1] = j;
                xi[2] = k;
                x.push_back(xi);
            }
        }
    }

    //set up a periodic box 10x10x10
    box_t box;
    box.period = vector<double>(3, 10.0);
    box.boxlo = vector<double>(3, 0.0);
    box.boxhi = vector<double>(3, 10.0);
    box.periodic = vector<bool>(3, true);
    
    //set up the interaction range (neighbor distance = 1.0)
    double tol = 0.01;
    double r0 = 1.0 + tol;

    //set up the default inputs for the RMS shape descriptor
    rmsdesc_info rms_info;

    //initialize a list to hold the clusters
    clstrlist_t clstrs;
    shpdesclist_t descriptors;
    
    //center the clusters at their centroid (default)
    localorder_info l_info;
    l_info.origin = CENTROID;

    localorder(x, box, r0, rmsdesc, &rms_info, clstrs, descriptors, l_info);
        
    //make a simple cubic reference cluster
    coordlist_t x_sc;
    for (int i=-1; i<=1; i++) {
        for (int j=-1; j<=1; j++) {
            for (int k=-1; k<=1; k++) {
                if (i*i + j*j + k*k < r0) {
                    coord_t xi(3);
                    xi[0] = i;
                    xi[1] = j;
                    xi[2] = k;                    
                    x_sc.push_back(xi);
                }
            }
        }
    }
    //create an RMS descriptor for the simple cubic reference structure
    shpdesc_t descriptor_sc;
    shapedata_t shape_sc(x_sc);
    rmsdesc(shape_sc, descriptor_sc, &rms_info);
    
    if (descriptors.size() != x.size()) {
        cerr << "local order create descriptors failed\n";
        errors++;
    }

    if (clstrs.size() != x.size()) {
        cerr << "local order create clusters failed\n";
        errors++;
    }

    for (unsigned int i=0; i<clstrs.size(); i++) {
        if (clstrs[i].size() != 7) {
            cerr << "local order clustering failed for cluster";
            cerr << " " << i << "\n";
            errors++;
        }
    }
        
    rmsassign_info assigninfo;
    for (unsigned int i=0; i<descriptors.size(); i++) {
        rmsfastassign(descriptors[i], descriptor_sc, &assigninfo);
        if (matchfundiff(descriptors[i], descriptor_sc, 0.0, 1.0) != 1.0) {
            cerr << "local order identify clusters failed for cluster";
            cerr << " " << i << "\n";
            errors++;
        }
    }
    
    if (errors) {
        cout << errors << " error(s) detected.\n";        
    }
    else {
        cout << "No errors detected.\n";
    }
    return 0;
}