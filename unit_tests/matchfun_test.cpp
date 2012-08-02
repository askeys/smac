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
 *  matchfun_test.cpp
 *  smac
 *
 *  Created by askeys on 9/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "smac.h"

using namespace std;
using namespace smac;

double tol = 1e-4;

bool different(double a, double b)
{
    if (fabs(a-b) > tol) {
        return true;
    }
    else {
        return false;
    }
}

void testall(vector<matchfun_t> matchfunlist, vector<string> nameofmatchfun, 
    int& errors, char* sd_type, shpdesc_t descriptor1, shpdesc_t descriptor2)
{
    for (unsigned int i=0; i<matchfunlist.size(); i++) {
        double match = matchfunlist[i](descriptor1, descriptor2, 0.0, 1.0);
        if (different(match, 1.0)) {
            cerr << "***ERROR: matching function " << nameofmatchfun[i];
            cerr << " failed test \"match " << sd_type <<" descriptors\"\n";
            cerr << "Ideal match = 1, ";
            cerr << "Actual match: " << match << endl;
            errors++;
        }
    }
    
    for (unsigned int i=0; i<matchfunlist.size(); i++) {
        double match = matchfunlist[i](descriptor1, descriptor2, 0.0, 100.0);
        if (different(match, 100.0)) {
            cerr << "***ERROR: matching function " << nameofmatchfun[i];
            cerr << " failed test \"match " << sd_type <<" descriptors\"";
            cerr << " with renormalization [0, 100]\n Ideal match = 100, ";
            cerr << "Actual match: " << match << endl;
            errors++;
        }
    }
    
    for (unsigned int i=0; i<matchfunlist.size(); i++) {
        double match = matchfunlist[i](descriptor1, descriptor2, 0.0, -1.0);
        if (different(match, -1.0)) {
            cerr << "***ERROR: matching function " << nameofmatchfun[i];
            cerr << " failed test \"match " << sd_type <<" descriptors\"";
            cerr << " with inverted renormalization [0, -1]\n Ideal match = -1\n";
            cerr << "Actual match: " << match << endl;
            errors++;
        }
    }
    
    //now negate one of the descriptors
    for (unsigned int i=0; i<descriptor2.size(); i++) {
        descriptor2[i] = -1.0*descriptor2[i];
    }
    
    for (unsigned int i=0; i<matchfunlist.size(); i++) {
        double match = matchfunlist[i](descriptor1, descriptor2, 0.0, 1.0);
        if (different(match, 0.0)) {
            cerr << "***ERROR: matching function " << nameofmatchfun[i];
            cerr << " failed test \"match negated " << sd_type <<" descriptors\"\n";
            cerr << "Ideal match = 0, ";
            cerr << "Actual match: " << match << endl;
            errors++;
        }
    }
    
    for (unsigned int i=0; i<matchfunlist.size(); i++) {
        double match = matchfunlist[i](descriptor1, descriptor2, -1.0, 100.0);
        if (different(match, -1.0)) {
            cerr << "***ERROR: matching function " << nameofmatchfun[i];
            cerr << " failed test \"match negated " << sd_type <<" descriptors\"";
            cerr << " with renormalization [-1, 100]\n Ideal match = -1, ";
            cerr << "Actual match: " << match << endl;
            errors++;
        }
    }
    
    for (unsigned int i=0; i<matchfunlist.size(); i++) {
        double match = matchfunlist[i](descriptor1, descriptor2, 1.0, -1.0);
        if (different(match, 1.0)) {
            cerr << "***ERROR: matching function " << nameofmatchfun[i];
            cerr << " failed test \"match negated " << sd_type <<" descriptors\"";
            cerr << " with inverted renormalization [1, -1]\n Ideal match = 1\n";
            cerr << "Actual match: " << match << endl;
            errors++;
        }
    }

}

int main(int argc, char** argv)
{
    cout << "Matching Function Test\n";
    cout << "---------------\n";

    int errors = 0;
    
    //make a list of all matching funcitons to loop over
    vector<matchfun_t> matchfunlist;
    vector<string> nameofmatchfun;
    matchfunlist.push_back(matchfundot);
    nameofmatchfun.push_back("matchfundot");
    matchfunlist.push_back(matchfundiff);
    nameofmatchfun.push_back("matchfundiff");
    matchfunlist.push_back(matchfundist);
    nameofmatchfun.push_back("matchfundist");
    matchfunlist.push_back(matchfunchisq);
    nameofmatchfun.push_back("matchfunchisq");

    //use nc components for each descriptor
    int nc = 20;
	shpdesc_t descriptor1(nc), descriptor2(nc);
    
    //try a real valued descriptor
    for (int i=0; i<nc; i++) {
        descriptor1[i] = component_t(drand48(), 0.0);
    }
    descriptor2 = descriptor1;
    
    testall(matchfunlist, nameofmatchfun, errors, "identical real-valued", descriptor1, descriptor2);

    for (int i=0; i<nc; i++) {
        descriptor1[i] = component_t(drand48(), drand48());
    }
    descriptor2 = descriptor1;

    testall(matchfunlist, nameofmatchfun, errors, "identical complex-valued", descriptor1, descriptor2);


    if (errors) {
        cout << errors << " error(s) detected.\n";        
    }
    else {
        cout << "No errors detected.\n";
    }
    return 0;

}

