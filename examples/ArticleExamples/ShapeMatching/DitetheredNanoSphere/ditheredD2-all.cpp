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
#include <utility>
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
		
	shpdistD2_info d2_args;
	d2_args.rmax = 4.0;    
	d2_args.nbin = 25;
    d2_args.box = box;
    
    shapedata_t s;
    xyzfile xyz1, xyz2, xyz3;
	shpdesc_t sd, sd_tet, sd_cyl, sd_dis;

    //we can see visually that this is most disordered here 
    //(and it obviously makes sense that it would be disordered at \epsilon ~ 0)
    //file = datapath + "74610000.dat";    
    loadxyz(filename4epsilon(0.1).c_str(), xyz1);
    x_dis = xyz1.x;
    vector<double> ones(xyz1.x.size(), 1.0);
    types_dis = xyz1.type;
	deltype(x_dis, types_dis, "C");
	deltype(x_dis, types_dis, "N");
    s.x = x_dis;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_dis.insert(sd_dis.end(), sd.begin(), sd.end());
    x_dis = xyz1.x;
    types_dis = xyz1.type;
	deltype(x_dis, types_dis, "O");
	deltype(x_dis, types_dis, "N");
    s.x = x_dis;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_dis.insert(sd_dis.end(), sd.begin(), sd.end());
    x_dis = xyz1.x;
    types_dis = xyz1.type;
	deltype(x_dis, types_dis, "O");
	deltype(x_dis, types_dis, "C");
    s.x = x_dis;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_dis.insert(sd_dis.end(), sd.begin(), sd.end());

    for (int i=0; i<sd_dis.size(); i++) {
        cout << sd_dis[i].real() << "\t";
    }
    cout << endl << endl;
    
    //we can see visually that this is a nice tetragonal phase
	loadxyz(filename4epsilon(0.8).c_str(), xyz2);
    x_tet = xyz2.x;
    types_tet = xyz2.type;
	deltype(x_tet, types_tet, "C");
	deltype(x_tet, types_tet, "N");
    s.x = x_tet;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_tet.insert(sd_tet.end(), sd.begin(), sd.end());
    x_tet = xyz2.x;
    types_tet = xyz2.type;
	deltype(x_tet, types_tet, "O");
	deltype(x_tet, types_tet, "N");
    s.x = x_tet;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_tet.insert(sd_tet.end(), sd.begin(), sd.end());
    x_tet = xyz2.x;
    types_tet = xyz2.type;
	deltype(x_tet, types_tet, "O");
	deltype(x_tet, types_tet, "C");
    s.x = x_tet;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_tet.insert(sd_tet.end(), sd.begin(), sd.end());

    for (int i=0; i<sd_tet.size(); i++) {
        cout << sd_tet[i].real() << "\t";
    }
    cout << endl << endl;

    //we can see visually that this is a nice cylinders phase
	loadxyz(filename4epsilon(1.6).c_str(), xyz3);
    x_cyl = xyz3.x;
    types_cyl = xyz3.type;
	deltype(x_cyl, types_cyl, "C");
	deltype(x_cyl, types_cyl, "N");
    s.x = x_cyl;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_cyl.insert(sd_cyl.end(), sd.begin(), sd.end());
    x_cyl = xyz3.x;
    types_cyl = xyz3.type;
	deltype(x_cyl, types_cyl, "O");
	deltype(x_cyl, types_cyl, "N");
    s.x = x_cyl;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_cyl.insert(sd_cyl.end(), sd.begin(), sd.end());
    x_cyl = xyz3.x;
    types_cyl = xyz3.type;
	deltype(x_cyl, types_cyl, "O");
	deltype(x_cyl, types_cyl, "C");
    s.x = x_cyl;
    s.f = ones;
	shpdistD2(s, sd, &d2_args);
    sd_cyl.insert(sd_cyl.end(), sd.begin(), sd.end());

    for (int i=0; i<sd_cyl.size(); i++) {
        cout << sd_cyl[i].real() << "\t";
    }
    cout << endl << endl;
            
    cout << "epsilon\tdis\thex\tcyl\n";
    double step = 0.05;
	for (double epsilon=0.1; epsilon<=1.61; epsilon+=step) {
        shpdesc_t sd_i;
        xyzfile xyz;
		loadxyz(filename4epsilon(epsilon).c_str(), xyz);
        coordlist_t x = xyz.x;
        vector<string> types =  xyz.type;
        x = xyz.x;
        types = xyz.type;
        deltype(x, types, "C");
        deltype(x, types, "N");
        s.x = x;
        s.f = ones;
        shpdistD2(s, sd, &d2_args);
        sd_i.insert(sd_i.end(), sd.begin(), sd.end());
        x = xyz.x;
        types = xyz.type;
        deltype(x, types, "O");
        deltype(x, types, "N");
        s.x = x;
        s.f = ones;
        shpdistD2(s, sd, &d2_args);
        sd_i.insert(sd_i.end(), sd.begin(), sd.end());
        x = xyz.x;
        types = xyz.type;
        deltype(x, types, "O");
        deltype(x, types, "C");
        s.x = x;
        s.f = ones;
        shpdistD2(s, sd, &d2_args);
        sd_i.insert(sd_i.end(), sd.begin(), sd.end());
		        
        //note: take a 3-pt average of the data to smooth it a little afterwards
        //this helps make the trends a little more clear
        //really, we should use more data, but to keep disk space down for the 
        //examples, we only use one dataset per epsilon
		cout << epsilon << "\t" << matchfundist(sd_i, sd_dis) << "\t"
			<< matchfundist(sd_i, sd_tet) << "\t"
			<< matchfundist(sd_i, sd_cyl) << "\n";
        
        if (epsilon >= 0.7) {
            step = 0.1;
        }
	}
	
	return 0;
}
