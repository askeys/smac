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

#define NINITALROTATIONS 100
#define NITER   10
#define QSTEP 0.025
#define BETA 20
#define TERMINATE -100.0

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

void flattencluster(coordlist_t& x)
{
    for (int i=0; i<x.size(); i++) {
        x[i][2] = 0.0;
    }
}

struct orientcluster_info
{
    orientcluster_info() : niter(30), nstartpoints(1000), dqmax(0.02), beta(20) {}
    
    int niter;
    int nstartpoints;
    double dqmax;
    double beta;
    fourier_info* fourierargs;
};

double computeangle(coordlist_t& x_orig, fourier_info& info)
{
    coordlist_t x = x_orig;
    flattencluster(x);
    clstrlist_t clusters;
    for (int j=0; j<x.size(); j++) {
        for (int k=0; k<3; k++) {
            x[j][k] -= x[0][k];
        }
    }
    shpdesc_t sd;
    shapedata_t shape(x);
    fourierdesc2d(shape, sd, &info);
    double xt = sd[0].real();
    double yt = sd[0].imag();
    double ltheta = atan2(yt, xt);
    return (-ltheta/(double)info.frequency[0]);

    /*
    double best_angle = 0.0;
    double min_energy = 1e100;
    double step = M_PI/100.0;
    vector<double> zaxis(3, 0.0);
    zaxis[2] = 1.0;
    for (int i=0; i<100; i++) {
        coordlist_t x = x_orig;
        rotateaa(x, i*step, zaxis);
        flattencluster(x);
        clstrlist_t clusters;
        for (int j=0; j<x.size(); j++) {
            for (int k=0; k<3; k++) {
                x[j][k] -= x[0][k];
            }
        }
        shpdesc_t sd;
        shapedata_t shape(x);
        fourierdesc2d(shape, sd, &info);
        double energy = 0.0;
        for (int k=0; k<sd.size(); k++) {
            energy -= sd[k].real();
        }
        
        if (energy < min_energy) {
            min_energy = energy;
            best_angle = i*step;
        }
    }
    cerr << best_angle << endl;
    rotateaa(x_orig, best_angle, zaxis);
    */
}

void orientcluster(coordlist_t& x_orig, vector<double>& q_best, orientcluster_info& info)
{
    q_best = vector<double>(4, 0.5);
    double e_best = 1e100;
    
    for (int s=0; s<info.nstartpoints; s++) {
        vector<double> q(4);
        randquaternion(q);
        double e_best_i = 0.0;
        vector<double> q_best_i = q;
        double dqmax = info.dqmax;
        for (int i=0; i<info.niter; i++) {
            vector<double> qr(4);
            q = q_best_i;
            randquaternion(qr);
            q[0] += dqmax*qr[0];
            q[1] += dqmax*qr[1];
            q[2] += dqmax*qr[2];
            q[3] += dqmax*qr[3];
            double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
            q[0] /= norm;
            q[1] /= norm;
            q[2] /= norm;
            q[3] /= norm;
            
            coordlist_t x = x_orig;
            rotateq(x, q);
            flattencluster(x);
            clstrlist_t clusters;
            for (int j=0; j<x.size(); j++) {
                for (int k=0; k<3; k++) {
                    x[j][k] -= x[0][k];
                }
            }
            shpdesc_t sd;
            shapedata_t shape(x);
            fourierdesc2d(shape, sd, info.fourierargs);
                
            //for (int k=0; k<sd.size(); k++) {
            //    cout << k << "\t" << abs(sd[k]) << endl;
            //}
                        
            double energy = 0.0;
            for (int k=0; k<sd.size(); k++) {
                energy -= abs(sd[k]);
            }
            //cout << i << "\t" << energy << "\t" << e_best_i << "\t" << exp(BETA*(e_best_i-energy)) << "\n";
            if (exp((info.beta+i)*(e_best_i-energy)) >= drand48()) {
                //cout << "accepted\n";
                e_best_i = energy;
                q_best_i = q;
            }
            else {
                dqmax *= 0.95;
                //cout << "rejected\n";

            }  
        }
        if (e_best_i < e_best) {
            e_best = e_best_i;
            q_best = q_best_i;
        }
        if (e_best_i < TERMINATE) {
            break;
        }
    }    
}

int main(int argc, char** argv)
{		
    system("rm tmp.xyz");
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/diffraction/structures/";
    string file = datapath + "sigma.xyz";
    
    double rnbrsq = 3.4*3.4;
    xyzfile xyz;
    loadxyz(file.c_str(), xyz);
    coordlist_t x = xyz.x;
    
    fourier_info info;
    info.shells.clear();
    info.shells.push_back(0.3);
    info.shells.push_back(1e100);

	//vector<int> k(1, 12);
    vector<int> k(13);
    for (int i=0; i<13; i++) {
        k[i] = i;
    }

    /*
    vector<int> k;
    k.push_back(3);
    k.push_back(4);
    k.push_back(6);
    k.push_back(8);
    k.push_back(12);
    */
    
	info.frequency = k;
    info.invariant = false;
        
    localorder_info lo_info;
    lo_info.origin = CENTER_PARTICLE;
    
    double L = 48.9997;
    box_t box;
    box.period = vector<double>(3, L);
    box.boxlo = vector<double>(3, -L/2.0);
    box.boxhi = vector<double>(3, L/2.0);

/*
    box.period = vector<double>(3, 20.0);
    box.boxlo = vector<double>(3, -10.0);
    box.boxhi = vector<double>(3, 10.0);
*/
    box.periodic = vector<bool>(3, false);
        
    coord_t origin(3, 0.0);
    
    vector<double> q(4, 0.0);
    //for (int i=0; i<1000; i++) { rand(); drand48();}
    
    clstrrulerange_info args;
    args.x = x;
    args.box = box;
    args.dim = 3;
    args.rcutsq = 1.5*1.5;
    clstrlist_t cluster;
    clstrshortrange(x, box, 1.5, clstrrulerange, &args, cluster);

    coordlist_t x_cluster;
    int ii;
    
    for (int i=0; i<500; i++) {
    //    double tmp = drand48();
      //  int tmpi = rand();
    }
    
    do {
        do {
            ii = rand()%x.size();
            //cerr << ii << endl;
        } while (cluster[ii].size() != 15);
    } while (distancesq(x[ii], origin, box.period, box.periodic) > 3.0*3.0);
    x_cluster.push_back(x[ii]);
    for (int jj=0; jj<x.size(); jj++) {
        if (distancesq(x[ii], x[jj], box.period, box.periodic) < rnbrsq) {
            //if (cluster[ii].size() == 15) {
                x_cluster.push_back(x[jj]);
            //}
        }
    }
    cout << x_cluster.size() << endl;
    cout << "orienting cluster\n";
    
    registersymmetryaxis_info reginfo;
    reginfo.niter = 30;
    reginfo.nstartpoints = 2000;
    reginfo.fourierargs = &info;

    registersymmetryaxis(x_cluster, q, reginfo);
    rotateq(x_cluster, q);
    rotateq(x, q);
    cout << "finding phase angle\n";
//    double angle = computeangle(x_cluster, info);
    vector<double> zaxis(3, 0.0);
    zaxis[2] = 1.0;
//    rotateaa(x, angle, zaxis); 

    vector<double> yaxis(3, 0.0);
    yaxis[1] = 1.0;

    //rotateaa(x, M_PI/2.0, yaxis);
    printxyz("tmp.xyz", x);
    
    cout << "printing diffraction image\n";
   // printxyz("tmp.xyz", x);
    printdiffractionpattern("matrix.txt", x);
    return 0;
}