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

// The point of this example is to demonstrate simple matching for a more 
// complex system.  The basic problem is we have an FCC crystal with HCP 
// stacking faults.  We use a rotation-invariant fourier descriptor to solve the
// problem, which makes the code somewhat shorter and more straightforward than 
// in the original stacking fault code

// Since we are only demonstrating how to match structures, we wont bother
// with registration, or anything fancy

#include <iostream>
#include <fstream>
#include <smac/smac.h>

using namespace std;
using namespace smac;

//a custom type to organize data for a reference structure
struct ref_structure
{
    coordlist_t x; 
    shpdesc_t descriptor;                   
    string name;                           
};

int main(int argc, char** argv)
{	    
    //load data from the example data directory
    //(note: you have to run the script in that directory to get the data 
    //from the server)
    
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/";        
    string file = datapath + "stacking/stackingfault.lammpstrj";
    lammpstrjfile lammpsdata;
    loadlammpstrj(file.c_str(), lammpsdata);
    coordlist_t x = lammpsdata.x;
    box_t box = lammpsdata.box;
    
    //make sure all the particles are in the box
    assert(inbox(x, box));
    
    //set up the interaction range
    double r0 = 1.3;

    //set up the default inputs for the RMS shape descriptor
    fourier_info finfo;
    finfo.frequency.clear();
    finfo.frequency.push_back(4);
    finfo.frequency.push_back(6);
    finfo.shells.clear();
    finfo.shells.push_back(0.1);
    finfo.shells.push_back(100);
    finfo.invariant = INVARIANT_Q | INVARIANT_W;
    
    //initialize a list to hold the clusters
    clstrlist_t clstrs;
    shpdesclist_t descriptors;
    
    //compute local clusters
    localorder_info l_info;
    //center the clusters at their centroid (default)
    l_info.origin = CENTROID;
    localorder(x, box, r0, fourierdesc3d, &finfo, clstrs, descriptors, l_info);

    //load data from the example data directory
    //(note: you have to run the script in that directory to get the data 
    //from the server)
    coordlist_t x_fcc, x_hcp;        
    file = datapath + "fcc_13.txt";
    load(file.c_str(), x_fcc);
    file = datapath + "hcp_13.txt";
    load(file.c_str(), x_hcp);
        
    shapedata_t shape_fcc(x_fcc), shape_hcp(x_hcp);
    shpdesc_t descriptor_fcc, descriptor_hcp;
    
    fourierdesc3d(shape_fcc, descriptor_fcc, &finfo);
    fourierdesc3d(shape_hcp, descriptor_hcp, &finfo);

    shpdesclist_t reference_library;
    reference_library.push_back(descriptor_hcp);
    reference_library.push_back(descriptor_fcc);
    
    map<int, string> name_of_reference_structure;
    name_of_reference_structure[0] = "HCP";
    name_of_reference_structure[1] = "FCC";
    
    matchlist_t match;    
    for (unsigned int i=0; i<descriptors.size(); i++) {
        bestmatch(descriptors[i], reference_library, matchfundiff, match);
        cout << i << "\t" << name_of_reference_structure[match[0].index] << "\t" << match[0].match << "\n";
    }

    //obviously, this could be in the loop above, but for clarity do it 
    //separately:
    ofstream f("stacking_color_by_structure.xyz");
    f << x.size() << "\n\n";
    for (unsigned int i=0; i<descriptors.size(); i++) {
        bestmatch(descriptors[i], reference_library, matchfundiff, match);
        if (match[0].index == 0) {
            f << "C\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n"; 
        }
        else {
            f << "N\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";         
        }
    }
    f.close();
    
    return 0;
}