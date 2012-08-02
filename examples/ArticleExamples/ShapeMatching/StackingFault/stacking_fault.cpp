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
// stacking faults.  We decide to use an RMS descriptor for consistency
// with part 1 of this example.  (Note however, that in general using rotation-
// invariant fourier descriptors might be a better idea for this type of 
// problem).  We can easily solve the problem with RMS descriptors, since we
// know there are only two types of clusters and since the structure is a 
// crystal, they must all be aligned.  First we match local structures within
// the crystal to label each atom type 1 or type 2.  Then we figure out what
// type 1 and type 2 stands for by matching with FCC and HCP reference clusters.
// This final step requires registration, but it is inexpensive for only two 
// structures.  In general the final step isn't really necessary; if we can
// label our structures type 1 and type 2, we can figure out what the types
// stand for via visual inspection.

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
    rmsdesc_info rms_info;

    //initialize a list to hold the clusters
    clstrlist_t clstrs;
    shpdesclist_t descriptors;
    
    //compute local clusters
    localorder_info l_info;
    //center the clusters at their centroid (default)
    l_info.origin = CENTROID;
    localorder(x, box, r0, rmsdesc, &rms_info, clstrs, descriptors, l_info);

    //compile a list of local structures on the fly 
    //add new structures to the list as we find them
    //determine if structures are new by matching with structures already in 
    //our list        
    vector<ref_structure> otf_reflib; 
    
    //list of ids for each reference cluster
    //id -> index in otf_reflib list
    vector<int> id_clstr(descriptors.size());                    

    rmsassign_info assigninfo;
    //if the match is better than this, it must be the same structure
    //to pick the value, we use some intuition from part one although trial and 
    //error would also work... 
    double min_match = 0.95;
    for (unsigned int i=0; i<descriptors.size(); i++) {
        bool found_match = false;
        for (unsigned int j=0; j<otf_reflib.size(); j++) {
            rmsfastassign(descriptors[i], otf_reflib[j].descriptor, &assigninfo);
            if (matchfundiff(descriptors[i], otf_reflib[j].descriptor) > min_match) {
                found_match = true;
                id_clstr[i] = j;
                break;
            }
        }
        if (!found_match) {
            id_clstr[i] = otf_reflib.size();
            coordlist_t x_clstr_i = subset(x, clstrs[i]);
            coord_t xc = x_clstr_i[0];
            unmap(x_clstr_i, xc, box);
            normusphere(x_clstr_i);
            ref_structure s;
            s.x = x_clstr_i;
            s.descriptor = descriptors[i];
            otf_reflib.push_back(s);            
        }        
    }
    int n_unknown = otf_reflib.size();
        
    //now we can figure out the true identify of the different clusters by 
    //matching our otf reference library with a true reference library 
    //we know our system is HCP / FCC so we load those two reference structures
    
    coordlist_t x_fcc, x_hcp;
    file = datapath + "fcc_13.txt";
    load(file.c_str(), x_fcc);
    file = datapath + "hcp_13.txt";
    load(file.c_str(), x_hcp);
    normusphere(x_fcc);
    normusphere(x_hcp);
    shapedata_t shape_fcc(x_fcc), shape_hcp(x_hcp);
    
    //compute descriptors for FCC and HCP reference clusters
    int n_ref = 2;
    vector<coordlist_t> x_ref(n_ref);
    vector<string> name_ref(n_ref);
    vector<shpdesc_t> descriptor_ref(n_ref);
    x_ref[0] = x_fcc;
    x_ref[1] = x_hcp;
    name_ref[0] = "FCC";
    name_ref[1] = "HCP";
    rmsdesc_info sdinfo;
    registericp_info register_info;
    register_info.maxiter = 10;
    register_info.maxrotations = 1024;
    register_info.tol = 0.01;
    
    //in V2 of this example, we use the bestmatch function, but here we do
    //matching manually.
    //there are only two unknowns and two reference structures, but let's write
    //it in a general way with loops for more re-usable code
    for (int i=0; i<otf_reflib.size(); i++) {
        min_match = 0.0;
        for (int j=0; j<n_ref; j++) {
            registericp(x_ref[j], otf_reflib[i].x, &register_info);
            shapedata_t shape_ref(x_ref[j]);
            shpdesc_t descriptor_ref;
            rmsdesc(shape_ref, descriptor_ref, &sdinfo);
            rmsfastassign(descriptor_ref, otf_reflib[i].descriptor, &assigninfo);
            double match = matchfundiff(descriptor_ref, otf_reflib[i].descriptor);

            if (match > min_match) {
                otf_reflib[i].name = name_ref[j];
                min_match = match;
            }
            if (min_match > 0.99) {
                break;
            }
        }
    }
    
    //print out the results
    for (unsigned int i=0; i<x.size(); i++) {
        cout << i << "\t" << otf_reflib[id_clstr[i]].name << "\n";
    }

    ofstream f("stacking_color_by_structure.xyz");
    f << x.size() << "\n\n";
    for (unsigned int i=0; i<x.size(); i++) {
        if (otf_reflib[id_clstr[i]].name == "FCC") {
            f << "C\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n"; 
        }
        else if (otf_reflib[id_clstr[i]].name == "HCP") {
            f << "N\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";         
        }
        else {
            f << "Cl\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";         
        }
    }
    f.close();

    return 0;
}