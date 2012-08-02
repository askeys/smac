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
 *  rmsdescriptor_test.h
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
    cout << "RMS Descriptor Test\n";
    cout << "---------------\n";

    int errors = 0;

    int np = 10;
    coord_t x11(3);
    coordlist_t x1(np,x11), x2(np,x11);
    for (int i=0; i<np; i++) {
        for (int j=0; j<3; j++) {
            x1[i][j] = drand48();
            x2[np-1-i][j] = x1[i][j];
        }
    }

    shapedata_t shape1(x1), shape2(x2);
    shpdesc_t descriptor1, descriptor2;
    rmsdesc_info sdinfo;

    rmsdesc(shape1, descriptor1, &sdinfo);
    rmsdesc(shape2, descriptor2, &sdinfo);

    double diff = 0.0;
    for (unsigned int ii=0; ii<descriptor1.size(); ii++) {
        int i=ii/3;
        int j=i%3;
        diff += (x1[i][j] - descriptor1[ii].real());
    }
    if (diff > 0.0) {
        cout << diff << "\n";
        cerr << "RMS basic test failed.\n";
        errors++;
    }

    if (matchfundiff(descriptor1, descriptor2, 0.0, 1.0) == 1.0) {
        cerr << "***WARNING: descriptors optimally correspond before ";
        cerr << "assignment.\n";
        errors++;
    }

    rmsassign_info assigninfo;
    rmsfastassign(descriptor1, descriptor2, &assigninfo);

    if (matchfundiff(descriptor1, descriptor2, 0.0, 1.0) != 1.0) {
        cerr << "RMS fastassign m=n failed\n";
        errors++;
    }

    //now try m>n
    rmsdesc(shape2, descriptor2, &sdinfo);
    descriptor2.push_back(component_t(drand48()));
    descriptor2.push_back(component_t(drand48()));
    descriptor2.push_back(component_t(drand48()));
    rmsfastassign(descriptor1, descriptor2, &assigninfo);
    descriptor2.resize(descriptor1.size());

    if (matchfundiff(descriptor1, descriptor2, 0.0, 1.0) != 1.0) {
        cerr << "RMS fastassign m>n failed\n";
        errors++;
    }

    //now try n>m
    rmsdesc(shape2, descriptor2, &sdinfo);
    descriptor2.push_back(component_t(drand48()));
    descriptor2.push_back(component_t(drand48()));
    descriptor2.push_back(component_t(drand48()));
    rmsfastassign(descriptor1, descriptor2, &assigninfo);
    descriptor2.resize(descriptor1.size());

    if (matchfundiff(descriptor2, descriptor1, 0.0, 1.0) != 1.0) {
        cerr << "RMS fastassign n>m failed\n";
        errors++;
    }

#ifdef _USE_HUNGARIAN
    rmsdesc(shape2, descriptor2, &sdinfo);
    rmsassign(descriptor1, descriptor2, &assigninfo);

    if (matchfundiff(descriptor1, descriptor2, 0.0, 1.0) != 1.0) {
        cerr << "RMS assign step failed\n";
        errors++;
    }
#endif

    if (errors) {
        cout << errors << " error(s) detected.\n";
    }
    else {
        cout << "No errors detected.\n";
    }
    return 0;
}
