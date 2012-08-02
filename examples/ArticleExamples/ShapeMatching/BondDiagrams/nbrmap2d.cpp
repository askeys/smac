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
#include <fstream>
#include <string>
#include <smgz/smgz.h>

#include <vmdstream/vmdstream.h>

using namespace std;
using namespace smgz;

void cubiclattice(coordlist_t& x)
{
    double p = 0.1;
    x.clear();
    for (double i=-5.0; i<5.0; i+=1.0) {
        for (double j=-5.0; j<5.0; j+=1.0) {
            for (double k=-5.0; k<5.0; k+=1.0) {
                coord_t xi(3);
                xi[0] = i + p*(drand48()-0.5);
                xi[1] = j + p*(drand48()-0.5);
                xi[2] = k + p*(drand48()-0.5);
                x.push_back(xi);
            }
        }
    }
}

void squarelattice(coordlist_t& x)
{
    double p = 0.4;
    x.clear();
    for (double i=-2.0; i<=2.0; i+=1.0) {
        for (double j=-2.0; j<=2.0; j+=1.0) {
            coord_t xi(3);
            xi[0] = i + p*(drand48()-0.5);
            xi[1] = j + p*(drand48()-0.5);
            xi[2] = 0.0;
            x.push_back(xi);
        }
    }
}

int main(int argc, char** argv)
{
    coordlist_t x;
    squarelattice(x);

    //set up a periodic box:
    box_t box;
/*
    box.period = vector<double>(3, 10.0);
    box.boxlo = vector<double>(3, -5.0);
    box.boxhi = vector<double>(3, 5.0);
*/
    double L = 5;
    box.period = vector<double>(3, L);
    box.boxlo = vector<double>(3, -L/2.);
    box.boxhi = vector<double>(3, L/2.);

    box.periodic = vector<bool>(3, true);

    for (int i=0; i<x.size(); i++) {
        pbc(x[i], box.period, box.periodic);
    }
    
    //make sure all the particles are in the box
    assert(inbox(x, box));
    
    //set up the interaction range
    double r0 = 1.5;

    shapedata_t input(x), output;
    neighbormap(input, 0.0, r0, box, output);
    normuball(output.x);
    
    int nsample = 200;
    ofstream f("bod_temp.xyz");
    f << nsample << "\n\n";
    for (int ii=0; ii<nsample; ii++) {
        int i = rand()%output.x.size();
        f << "O\t" << output.x[i][0] << "\t" << output.x[i][1] << "\t" << output.x[i][2] << "\n";
    }
    f.close();
    
	vmdsock_t vmdsock = newvmdsock();
	vmdstream vmd(vmdsock); 
    
    vmd << "molecule load xyz bod_temp.xyz\n";
    vmd << "set sel [atomselect top all]" << "\n";
    
    vmd << "set molid 0\n";
    vmd << "mol delrep 0 $molid\n";
    vmd << "mol addrep $molid\n";

    vmd << "menu main move 700 100\n";
    vmd << "display projection orthographic\n";
    vmd << "axes location off\n";
    vmd << "mol modstyle 0 0 Points 1.000000\n";
    vmd << "mol modcolor 0 0 ColorID 16\n";

    //vmd << "mol modcolor 0 0 Chain\n";

    vmd << "mol new\n";		
    vmd << "graphics 1 color white\n";		
    vmd << "graphics 1 materials on\n";		
    vmd << "graphics 1 material AOEdgy\n";		
    //cylinder {x1 y1 z1} {x2 y2 z2} [radius r] [resolution n] [filled $<$yes|no$>$]: 
    vmd << "graphics 1 cylinder {0 0 -0.2} {0 0 -0.1} radius 1.2 resolution 20 filled yes\n";
    vmd << "color Display Background white\n";
    
    vmd.flush();
    
	closevmdsock(vmdsock);
    
	return 0;
}
