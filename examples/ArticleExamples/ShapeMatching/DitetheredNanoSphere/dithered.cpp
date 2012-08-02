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
 *  dithered.cpp
 *  glotzilla
 *
 *  Created by askeys on 9/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <smac/smac.h>

using namespace std;
using namespace smac;

// The purpose of this example is to show how a ditethered nano sphere system
// goes through two transition as a function of the potential attraction 
// \epsilon.  The order parameter that we use is based on the the distribution
// of local density maps projected on the surface of spheres at different 
// lengthscales throughout the sample.  The density maps are indexed as rotation 
// invariant fourier descriptors for matching.

// get the filename corresponding to a particular epsilon value
// (has to do with idiosyncracies in the directory structure for this
// particular dataset)
string filename4epsilon(double epsilon)
{
    char eps[8];
    sprintf(eps, "%1.3f", epsilon);
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/dtns/";
	string file;
    file = datapath + "epsilon=" + eps + "/final/*.dat";
    ostringstream lscommand;
    lscommand << "ls " << file << " > tmp.txt";
    system(lscommand.str().c_str());
    ifstream f("tmp.txt");
    f >> file;
    if (f.fail()) {
        cerr << "ERROR opening " << file << endl;
        exit(1); 
    }
    f.close();
    return file;
}

int main(int argc, char** argv)
{
	box_t box;	
	coordlist_t x_dis, x_tet, x_cyl;
	vector<string> types_dis, types_tet, types_cyl;

	box.period = coord_t(3, 26.111928);
	box.boxhi = coord_t(3, 26.111928/2.0);
	box.boxlo = coord_t(3, -26.111928/2.0);
	box.periodic = vector<bool>(3, true);	
		
	fourier_info f_args;
	f_args.invariant = INVARIANT_Q;
	f_args.frequency.clear();
    for (int i=4; i<=12; i++) {
        f_args.frequency.push_back(i);
    }
	f_args.shells.clear();
	f_args.shells.push_back(0.1);
	f_args.shells.push_back(2.8);
	f_args.shells.push_back(4.15);
	f_args.shells.push_back(5.4);
	double max_range = f_args.shells[f_args.shells.size()-1];
	
	shpdist_info d_args;
	d_args.nbin = 10;
	d_args.min = 0.0;
	d_args.max = 1.0;

    xyzfile xyz1, xyz2, xyz3;

    //we can see visually that this is most disordered here 
    //(and it obviously makes sense that it would be disordered at \epsilon ~ 0)
    //file = datapath + "74610000.dat";    
    loadxyz(filename4epsilon(0.1).c_str(), xyz1);
    x_dis = xyz1.x;
    types_dis = xyz1.type;
	deltype(x_dis, types_dis, "C");
	deltype(x_dis, types_dis, "N");
	
    //we can see visually that this is a nice tetragonal phase
	loadxyz(filename4epsilon(0.8).c_str(), xyz2);
    x_tet = xyz2.x;
    types_tet = xyz2.type;
	deltype(x_tet, types_tet, "C");
	deltype(x_tet, types_tet, "N");

    //we can see visually that this is a nice cylinders phase
	loadxyz(filename4epsilon(1.6).c_str(), xyz3);
    x_cyl = xyz3.x;
    types_cyl = xyz3.type;
	deltype(x_cyl, types_cyl, "C");
	deltype(x_cyl, types_cyl, "N");
    
    //now lets quantify what we see visually:
    
	shpdesc_t sd_tet, sd_cyl, sd_dis;
    shpdesclist_t sd_all;

    localorder_info lo_info;
    lo_info.origin = CENTROID;

    clstrlist_t clusters;
	localorder(x_dis, box, max_range, fourierdesc3d, &f_args, clusters, sd_all, lo_info);
	shpdist(sd_all, sd_dis, &d_args);
	//shpavg(sd_all, sd_dis);
		
	localorder(x_tet, box, max_range, fourierdesc3d, &f_args, clusters, sd_all, lo_info);	
	shpdist(sd_all, sd_tet, &d_args);
	//shpavg(sd_all, sd_tet);
	
	localorder(x_cyl, box, max_range, fourierdesc3d, &f_args, clusters, sd_all, lo_info);
	shpdist(sd_all, sd_cyl, &d_args);
	//shpavg(sd_all, sd_cyl);
	
    double step = 0.05;
	for (double epsilon=0.1; epsilon<=1.61; epsilon+=step) {
        xyzfile xyz;
		loadxyz(filename4epsilon(epsilon).c_str(), xyz);
        coordlist_t x = xyz.x;
        vector<string> types =  xyz.type;
		deltype(x, types, "C");
		deltype(x, types, "N");
		
		shpdesc_t sd;								
		localorder(x, box, max_range, fourierdesc3d, &f_args, clusters, sd_all, lo_info);
		shpdist(sd_all, sd, &d_args);
		//shpavg(sd_all, sd);
        
        //note: take a 3-pt average of the data to smooth it a little afterwards
        //this helps make the trends a little more clear
		cout << epsilon << "\t" << matchfundist(sd, sd_dis) << "\t"
			<< matchfundist(sd, sd_tet) << "\t"
			<< matchfundist(sd, sd_cyl) << "\n";
        
        if (epsilon >= 0.7) {
            step = 0.1;
        }
	}
	
	return 0;
}
