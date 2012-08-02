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

#define NINITALROTATIONS 300
#define NITER   20
#define QSTEP 0.01
#define BETA 20
#define TERMINATE -10000000.0

void printxyz(const char* filename, coordlist_t& x)
{
    ofstream f(filename, ios::app);
    f << x.size() << "\n\n";
    for (unsigned int i=0; i<x.size(); i++) {
        f << "H\t" << x[i][0] << "\t"<< x[i][1] << "\t"<< x[i][2] << "\n";
    }
    f.close();
}

void printdiffractionpattern(const char* filename, coordlist_t& x)
{
    double qmax = 15;
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

void flatten_system(coordlist_t& x, coordlist_t& x_flat)
{
    x_flat = x;
    for (int i=0; i<x_flat.size(); i++) {
        x_flat[i][2] = 0.0;
    }

    vector<bool> erased(x.size(), false);
    for (int i=0; i<x.size(); i++) {
        if (!erased[i]) {
            for (int j=i+1; j<x.size(); j++) {
                if (!erased[j]) {
                    if (distancesq(x_flat[i], x_flat[j]) < 0.2*0.2) {
                        erased[j] = true;
                    }
                }
            }
        }
    }
    x_flat.clear();
    for (int i=0; i<x.size(); i++) {
        if (!erased[i]) {
            coord_t xi = x[i];
            xi[2] = 0.0;
            x_flat.push_back(xi);
        }
    }
}

double cell_energy(coordlist_t& x, box_t& box, double cell_size, vector< vector<int> >& grid)
{
    double invbox0 = 1.0/box.period[0];
    double invbox1 = 1.0/box.period[1];
    
    int ncellx = box.period[0]/cell_size;
    int ncelly = box.period[1]/cell_size;
    
    for (int i=0; i<grid.size(); i++) {
        for (int j=0; j<grid[i].size(); j++) {
            grid[i][j] = 0;
        }
    }
    
    double energy = 0.0;
    for (int i=0; i<x.size(); i++) {
        int ii = ncellx*((x[i][0]-box.boxlo[0])*invbox0);
        int jj = ncelly*((x[i][1]-box.boxlo[1])*invbox1);
        if (!grid[ii][jj]) {
            grid[ii][jj] = 1;
            energy+=1.0;
        }
    }
    /*
    double energy = 0.0;
    for (int i=0; i<grid.size(); i++) {
        for (int j=0; j<grid[i].size(); j++) {
            for (int ii=-1; ii<=1; ii++) {
                for (int jj=-1; jj<=1; jj++) {
                    if (ii == 0 && jj==0) {
                        continue;
                    }
                    if (i == 0 || i == (ncellx-1) || j == 0 || j == (ncelly-1)) {
                        continue;
                    }
                    int delta = grid[i][j]-grid[i+ii][j+jj];
                    energy -= delta*delta;
                }
            }
        }
    }
    */
    
    return energy;
}

void find_face(coordlist_t& x_orig, box_t& box, vector<double>& qbest)
{    
    coordlist_t x = x_orig;
    double ebest = 1e100;
    double cell_size = 0.1;
    int ncellx = box.period[0]/cell_size;
    int ncelly = box.period[1]/cell_size;
    vector< vector<int> > grid(ncellx, vector<int>(ncelly, 0));
    
    for (int s=0; s<NINITALROTATIONS; s++) {
        vector<double> q(4);
        randquaternion(q);
        double ebest_i = 1e100;
        vector<double> qbest_i = q;
        for (int i=0; i<NITER; i++) {
            vector<double> qr(4);
            randquaternion(qr);
            q[0] += QSTEP*qr[0];
            q[1] += QSTEP*qr[1];
            q[2] += QSTEP*qr[2];
            q[3] += QSTEP*qr[3];
            double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
            q[0] /= norm;
            q[1] /= norm;
            q[2] /= norm;
            q[3] /= norm;
            
            x = x_orig;
            rotateq(x, q);
            
            double energy = cell_energy(x, box, cell_size, grid);
            cout << energy << "\t" << ebest_i << "\t" << exp(BETA*(ebest_i-energy)) << "\n";

            if (energy < ebest_i) {
                cout << "accepted\n";
                ebest_i = energy;
                qbest_i = q;
            }
            else {
                q = qbest_i;
                cout << "rejected\n";
            }
        }
        if (ebest_i < ebest) {
            qbest = qbest_i;
            ebest = ebest_i;
        }
        if (ebest_i < TERMINATE) {
            break;
        }
    }    
}

void planedescriptor(coordlist_t& x, box_t& box, double r_nbr, shpdesc_t& sd_avg)
{
    fourier_info info;
	vector<int> k(13);
	for (int i=0; i<13; i++) {
		k[i] = i;
	}
	info.frequency = k;
    info.invariant = false;
    
    localorder_info lo_info;
    lo_info.origin = CENTER_PARTICLE;

    shpdesclist_t sd_all;
    clstrlist_t clusters;
    localorder(x, box, r_nbr, fourierdesc2d, &info, clusters, sd_all, lo_info);

    for (int i=0; i<clusters.size(); i++) {
        cout << clusters[i].size() << endl;
    }

    sd_avg.resize(sd_all[0].size(), 0.0);
    for (int j=0; j<sd_all.size(); j++) {
        for (int k=0; k<sd_avg.size(); k++) {
            sd_avg[k] += sd_all[j][k];
        }
    }
    for (int k=0; k<sd_avg.size(); k++) {
        sd_avg[k] /= (double)sd_all.size();
    }
}


int main(int argc, char** argv)
{		
    system("rm tmp*.xyz");
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/fcc_crystal/";
    string file = datapath + "dump.xyz";
    
    double r_nbr = 1.2;
    xyzfile xyz;
    loadxyz(file.c_str(), xyz);
    coordlist_t x = xyz.x;
        
    box_t box;
    box.period = vector<double>(3, 20.0);
    box.boxlo = vector<double>(3, -10.0);
    box.boxhi = vector<double>(3, 10.0);
    box.periodic = vector<bool>(3, false);

    coordlist_t x_orig = x;
    
    vector<double> q_face(4);
    find_face(x, box, q_face);
    
    //now that we have a face, align the face
    shpdesc_t descriptor;
          /*
        int k = max_symmetry(descriptors);
        double real = descriptors[k];
        rotateaa(x, angle, zaxis);
    }
    */
    
    
    
    coordlist_t x_flat;
    x = x_orig;
    rotateq(x, q_face);
    flatten_system(x, x_flat);
    planedescriptor(x_flat, box, r_nbr, descriptor);
    cout << descriptor << endl;

    printxyz("tmp1.xyz", x);
    exit(1);
    x = x_orig;
    vector<double> yaxis(3, 0.0);
    yaxis[1] = 1.0;
    rotateq(x, q_face);
    //rotateaa(x, M_PI/2.0, yaxis);
    rotated(x, 0, 90, 0);
    printxyz("tmp2.xyz", x);
    
    x = x_orig;
    vector<double> zaxis(3, 0.0);
    zaxis[2] = 1.0;
    rotateq(x, q_face);
    rotated(x, 0, 0, 90);
    printxyz("tmp3.xyz", x);
    
    printdiffractionpattern("matrix.txt", x);
    return 0;
}