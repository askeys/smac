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
 *  test_sd.cpp
 *  glotzilla
 *
 *  Created by askeys on 9/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 *  To compile: g++ test_sd.cpp -Iglotzilla_path glotz_lib_path/libglotzilla.a 
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <smgz/smgz.h>
#include <vmdstream/vmdstream.h>

using namespace std;
using namespace smgz;

#define R0 0.5
#define R1 0.5

#define BONDCUT 0.992
#define FBONDCUT 0.5

void load2d(const char* filename, coordlist_t& x)
{
    x.clear();
    ifstream f(filename);
    assert(f.good());
    while (1) {
        double i, j;
        f >> i >> j;
        if (f.fail()) {
            break;
        }
        coord_t xi(3, 0.0);
        xi[0] = i;
        xi[1] = j;
        x.push_back(xi);
    }    
}

void load2d(const char* filename, coordlist_t& x, vector<int>& types)
{
    load2d(filename, x);
    types.resize(x.size());
    
    for (int i=0; i<x.size(); i++) {
        types[i] = 0;
        double minr = 100;
        for (int j=0; j<x.size(); j++) {
            if (i==j) {
                continue;
            }
            minr = min(distancesq(x[i], x[j]), minr);
        }
        if (minr > .11) {
            types[i] = 1;
        }
    }
}

int main(int argc, char** argv)
{		
    string datapath = "/Users/askeys/Data/GlotzerChandler/granular/granular_data_lg/N1904_AF773_U40/";
    string file = datapath + "10000.txt";
    
    coordlist_t x;
    vector<int> types;
    load2d(file.c_str(), x, types);
    transcentroid(x);
    
    double box_size = 10.0;
	box_t box;
	box.period = vector<double>(3, 2*box_size);
	box.boxhi = vector<double>(3, box_size);
	box.boxlo = vector<double>(3, -box_size);
	box.boxhi[2] = 3.0;
	box.boxlo[2] = -3.0;
	box.period[2] = 6;
	box.periodic[2] = false;

    assert(inbox(x, box));

	//set up the shape descriptor:
    //fourier descriptor, k=5,6,7,8,9,10, single radial shell
	rmsdesc_info sd_args;
	
    //set up parameters for the crystal analysis function
	crystalanalysis_info c_info;
	c_info.rangedescriptor = R0;
	c_info.rangecompare = R1;
	c_info.rangecluster = R0;
	c_info.descriptor = rmsdesc;
	c_info.descriptorarg = &sd_args;
	c_info.matchfun = matchfundotrms;
	c_info.crystalbondcut = BONDCUT;
	c_info.crystalparticlecut = FBONDCUT;
	
    //find the crystal grains:
    clstrlist_t crystal_grains;
    crystalgrains(x, box, c_info, crystal_grains);

    //open a socket to vmd (required vmdstream header file)
	vmdsock_t vmdsock = newvmdsock();
	vmdstream f(vmdsock); 
    int color = 0;
    vector<bool> in_crystal(x.size(), false);
    for (unsigned int i=0; i<crystal_grains.size(); i++) {
        f << "draw color " << color << endl;
        if (crystal_grains[i].size() > 3) {
            for (unsigned int j=0; j<crystal_grains[i].size(); j++) {
                f << "draw sphere {";
                f << x[crystal_grains[i][j]][0] << " ";
                f << x[crystal_grains[i][j]][1] << " ";
                f << x[crystal_grains[i][j]][2] << "} ";
                if (types[crystal_grains[i][j]] == 0) {
                    f << "radius 0.154" << endl;
                }
                else {
                    f << "radius 0.165" << endl;
                }
                in_crystal[crystal_grains[i][j]] = true;
            }
            color++;
        }
    }
    f << "draw color " << color << endl;
    for (unsigned int i=0; i<in_crystal.size(); i++) {
        if (!in_crystal[i]) {
            f << "draw sphere {" << x[i][0] << " " << x[i][1] << " " << x[i][2] << "} ";
            if (types[i] == 0) {
                f << "radius 0.154" << endl;
            }
            else {
                f << "radius 0.165" << endl;
            }
        }
    }
    f.flush();
	closevmdsock(vmdsock);	
}