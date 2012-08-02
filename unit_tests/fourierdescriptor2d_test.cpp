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
    cout << "Fourier Descriptor 2D Test\n";
    cout << "---------------\n";

    int errors = 0;

    vector<shapedata_t> poly;
    for (int i=0; i<11; i++) {
        coordlist_t x;
        polygon(i, x);
        shapedata_t s(x);
        poly.push_back(s);
    }
    
	fourier_info args;
	args.frequency = vector<int>(11);
    for (int i=0; i<11; i++) {
        args.frequency[i] = i;
    }
    
	args.invariant = INVARIANT_Q;

    for (int i=3; i<11; i++) {      
        shpdesc_t c;
        fourierdesc2d(poly[i], c, &args);
        if (!similar(c[0].real(), 1.0, 0.01)) {
            cerr << "***ERROR: bad value for cl polygon size " << i << " l = " << 0 << ": " << c[0] << "\n";
            errors++;    
        }
        for (int j=1; j<i; j++) {
            if (!similar(c[j].real(), 0.0, 0.01)) {
                cerr << "***ERROR: bad value for cl polygon size " << i << " l = " << j << ": " << c[j] << "\n";
                errors++;    
            }        
        }
        if (!similar(c[i].real(), 1.0, 0.01)) {
            cerr << "***ERROR: bad value for cl polygon size " << i << " l = " << i << ": " << c[i] << "\n";
            errors++;    
        }
    }
    
    //test tabulated trig functions
    
    args.trigtablesize = 1000;
    for (int i=3; i<11; i++) {      
        shpdesc_t c;
        fourierdesc2d(poly[i], c, &args);
        if (!similar(c[0].real(), 1.0, 0.01)) {
            cerr << "***ERROR: bad value for cl polygon size " << i << " l = " << 0 << ": " << c[0] << "\n";
            errors++;    
        }
        for (int j=1; j<i; j++) {
            if (!similar(c[j].real(), 0.0, 0.01)) {
                cerr << "***ERROR: bad value for cl polygon size " << i << " l = " << j << ": " << c[j] << "\n";
                errors++;    
            }        
        }
        if (!similar(c[i].real(), 1.0, 0.01)) {
            cerr << "***ERROR: bad value for cl polygon size " << i << " l = " << i << ": " << c[i] << "\n";
            errors++;    
        }
    }
    
    //test reconstruction
    /*
    args.nsector=12;
    shpdesc_t c6(1);
    args._frequency = 6;
    fouriercoeff2d(poly[6], c6, &args);
    
    weight_t f(12, 0.0);
    coordlist_t xtmp(12);
    shapedata_t shape(xtmp, f);
    
    invfouriercoeff2d(c6, shape, &args);
    
    for (int i=0; i<shape.f.size(); i++) {
        cout << shape.f[i] << endl;
    }
    */
    
    if (errors) {
        cout << errors << " error(s) detected.\n";        
    }
    else {
        cout << "No errors detected.\n";
    }
    return 0;
}