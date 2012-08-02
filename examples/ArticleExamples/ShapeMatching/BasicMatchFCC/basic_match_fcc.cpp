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
 *  basic_match_fcc.h
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

// The point of this example is to demonstrate simple matching
// Since we are only demonstrating how to match structures, we wont bother
// with registration, or anything fancy

#include <iostream>
#include <fstream>
#include <smac/smac.h>

using namespace std;
using namespace smac;

int main(int argc, char** argv)
{	
    double perturbation_size = 0.1;
    if (argc > 1) {
        perturbation_size = atof(argv[1]);
    }
    
    //load data from the example data directory
    //(note: you have to run the script in that directory to get the data 
    //from the server)
    coordlist_t x_unknown, x_fcc, x_ico, x_hcp;
    string file;
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/";
        
    file = datapath + "fcc_13.txt";
    load(file.c_str(), x_fcc);
    file = datapath + "hcp_13.txt";
    load(file.c_str(), x_hcp);
    file = datapath + "ico_13.txt";
    load(file.c_str(), x_ico);
        
    //set unknown = to fcc, except with a small perturbation    
    x_unknown = x_fcc;    
    for (int i=0; i<13; i++) {
        for (int j=0; j<3; j++) {
            x_unknown[i][j] += 2.0*perturbation_size*(drand48()-0.5);
        }
    }
    
    //register the coordinates beforehand
    //(note, if you use the optional arguments, you can make this step much faster)
    registericp(x_fcc, x_unknown);
    registericp(x_hcp, x_unknown);
    registericp(x_ico, x_unknown);

    shapedata_t shape_fcc(x_fcc), shape_hcp(x_hcp), shape_ico(x_ico), 
        shape_unknown(x_unknown);
    
    shpdesc_t descriptor_fcc, descriptor_hcp, descriptor_ico, 
        descriptor_unknown;
    
    rmsdesc_info sdinfo;
    rmsdesc(shape_fcc, descriptor_fcc, &sdinfo);
    rmsdesc(shape_hcp, descriptor_hcp, &sdinfo);
    rmsdesc(shape_ico, descriptor_ico, &sdinfo);
    rmsdesc(shape_unknown, descriptor_unknown, &sdinfo);

    shpdesclist_t reference_library;
    reference_library.push_back(descriptor_hcp);
    reference_library.push_back(descriptor_fcc);
    reference_library.push_back(descriptor_ico);
    
    map<int, string> name_of_reference_structure;
    name_of_reference_structure[0] = "HCP";
    name_of_reference_structure[1] = "FCC";
    name_of_reference_structure[2] = "ICO";
    
    rmsfastassign_info info;
    matchlist_t match;

    bestmatch(descriptor_unknown, reference_library, matchfundiff, 
        match, rmsfastassign, &info);
    
    cout << "Match values, best->worst:\n";
    for (int i=0; i<match.size(); i++) {
        cout << i+1 << " " << name_of_reference_structure[match[i].index];
        cout << "\t" << match[i].match << "\n";
    }
    
    return 0;
}
