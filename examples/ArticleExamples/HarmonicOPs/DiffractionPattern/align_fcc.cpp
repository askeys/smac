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
*/
/*
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
#include <smac/smac.h>

using namespace std;
using namespace smac;

void printxyz(const char* filename, coordlist_t& x)
{
    ofstream f(filename);
    f << x.size() << "\n\n";
    for (unsigned int i=0; i<x.size(); i++) {
        f << "H\t" << x[i][0] << "\t"<< x[i][1] << "\t"<< x[i][2] << "\n";
    }
    f.close();
}

void printdiffractionpattern(const char* filename, coordlist_t& x)
{
    double qmax = 30;
    double qstep = 0.25;
    int nq = 2*qmax/qstep+1;
    std::vector< std::vector< std::complex<double> > > complex_image(nq, std::vector< std::complex<double> >(nq, 0.0));
            
    for (unsigned int i=0; i<x.size(); i++) {
        for (int ii=0; ii<nq; ii++) {
            double qi = (ii-nq/2)*qstep;
            for (int jj=0; jj<nq; jj++) {
                double qj = (jj-nq/2)*qstep;
                complex_image[ii][jj] += std::complex<double>(cos(qi*x[i][0]+qj*x[i][1]), -sin(qi*x[i][0]+qj*x[i][1]));
            }
        }
    }        
    ofstream f(filename);
    for (unsigned int i=0; i<complex_image.size(); i++) {
        for (unsigned int j=0; j<complex_image[i].size(); j++) {
            f << abs(complex_image[i][j]) << "\t";
        }
        f << "\n";
    }
    f.close();
}

int main(int argc, char** argv)
{		
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/diffraction/structures/";
    string file = datapath + "fcc.xyz";
    
    xyzfile xyz;
    loadxyz(file.c_str(), xyz);
    coordlist_t x = xyz.x;

    for (int i=0; i<100; i++) {
        drand48();
    }
    
    //set a randome initial rotation:
    vector<double> q(4, 0.0);
    randquaternion(q);
    rotateq(x, q);
    //perturb the positions a little from the ideal ones:
    for (int i=0; i<x.size(); i++) {
        for (int k=0; k<3; k++) {
            x[i][k] += 0.1*randgauss();
        }
    }

    for (int i=0; i<200; i++) {
        drand48();
    }

    printxyz("fcc_original.xyz", x);  

    fourier_info info;
    info.shells.clear();
    info.shells.push_back(0.3);
    info.shells.push_back(1e100);
	info.frequency = vector<int>(1, 8);
	//info.frequency = vector<int>(1, 12);
    info.invariant = false;
            
    double L = 20.0;
    box_t box;
    box.period = vector<double>(3, L);
    box.boxlo = vector<double>(3, -L/2.0);
    box.boxhi = vector<double>(3, L/2.0);
    box.periodic = vector<bool>(3, false);
        
    coord_t origin(3, 0.0);    
    
    coordlist_t x_cluster;
    for (int i=0; i<x.size(); i++) {
        if (distancesq(x[i], origin, box.period, box.periodic) < 3.0*3.0) {
            x_cluster.push_back(x[i]);
        }
    }
    
    cout << "particles in test cluster: " << x_cluster.size() << endl;
    cout << "orienting cluster\n";
    
    registersymmetryaxis_info reginfo;
    reginfo.niter = 30;
    reginfo.nstartpoints = 1000;
    reginfo.fourierargs = &info;

    registersymmetryaxis(x_cluster, q, reginfo);
    cout << "finding phase angle\n";
    rotateq(x_cluster, q);
    double angle = registerzerocomplexcomonent(x_cluster, info);
    cout << "angle: " << angle << endl;
    vector<double> zaxis(3, 0.0);
    zaxis[2] = 1.0;
    rotateq(x, q);
    rotateaa(x, M_PI/4.0, zaxis);
    printxyz("fcc_register_step1.xyz", x);  
    rotateaa(x, angle, zaxis);
    printxyz("fcc_register_step2.xyz", x);  
    cout << "printing diffraction image\n";
    printdiffractionpattern("matrix.txt", x);
    return 0;
}