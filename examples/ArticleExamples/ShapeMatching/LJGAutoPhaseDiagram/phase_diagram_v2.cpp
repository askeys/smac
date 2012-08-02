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
#include <smac/smac.h>

using namespace std;
using namespace smac;

//a function to load the custom data format:
void load_special(string datafile, coordlist_t& x)
{
    ifstream is(datafile.c_str());
    if (is.fail()) {
        cerr << "Failed to open file " << datafile << endl;
        exit(1);
    }
    while (1) {
        string str;
        getline(is, str);
        if (str == "//data") {
            break;
        }
    }
    x.resize(1024, vector<double>(3, 0.0));
    for (int i=0; i<1024; i++) {
        string tmp;
        double xi, yi, zi;
        is >> tmp >> xi >> yi >> zi;
        x[i][0] = xi;
        x[i][1] = yi;
        x[i][2] = zi;
    }
    is.close();
}

int main(int argc, char** argv)
{		
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/ljg/";

    fourier_info info;
	vector<int> k(6);
	for (int i=0; i<=5; i++) {
		k[i] = i+5;
	}
	info.frequency = k;
    info.invariant = true;
        
    shpdesclist_t sd_all;
    localorder_info lo_info;
    lo_info.origin = CENTER_PARTICLE;
    
    box_t box;
    box.period = vector<double>(3, 60);
    box.boxlo = vector<double>(3, -30);
    box.boxlo = vector<double>(3, 30);
    box.periodic = vector<bool>(3, true);
    box.periodic[2] = false;
    box.boxlo[2] = -4;
    box.boxhi[2] = 4;
    box.period[2] = 8;

    //1) make a list of global descriptors for all statepoints

    map<string, map<string, shpdesc_t> > descriptor;
    //Note: use strings for the map key rather than doubles, since map<double> 
    //can act funny because == operator is not as well defined as for int, str    
    double r0_step = 0.01; //0.01;
    double eps_step = 0.1;
    double r_nbr = 1.5;
    
    for (double r0_i = 1.11; r0_i<=2.10; r0_i+=r0_step) {
        for (double eps_i = 0.1; eps_i<=5.0; eps_i+=eps_step) {
            ostringstream filename_i;
            char r0str[8], epsstr[8];
            sprintf(r0str, "%1.2f", r0_i);
            sprintf(epsstr, "%1.1f", eps_i);
            filename_i << datapath << "ljg_" << r0str << "_25_" << epsstr << ".pos";
            cout << filename_i.str() << endl;
            coordlist_t x_i; 
            load_special(filename_i.str(), x_i);
            clstrlist_t clusters;
            localorder(x_i, box, r_nbr, fourierdesc2d, &info, clusters, sd_all, lo_info);
            shpdesc_t sd_distribution;  //distribution of local descriptors
            shpdist_info dist_info;
            dist_info.nbin = 15;
            shpdist(sd_all, sd_distribution, &dist_info);
            descriptor[r0str][epsstr] = sd_distribution;
        }
    }
    
    //2) Match statepoints with their neighbors to determine which statepoints
    //   represent stable/changing regions of phase space 
    int nbr_range = 1;  //range for checking phase diagram neighbors
    ofstream matrix("matrix.txt");
    for (double eps_i = 0.1; eps_i<=5.0; eps_i+=eps_step) {
        for (double r0_i = 1.11; r0_i<=2.10; r0_i+=r0_step) {
            //loop over phase diagram neighbors
            double total = 0.0;
            int count = 0;
            for (double r0_j = r0_i-nbr_range*r0_step; r0_j<=r0_i+r0_step; r0_j+=nbr_range*r0_step) {
                for (double eps_j = eps_i-nbr_range*eps_step; eps_j<=eps_i+eps_step; eps_j+=nbr_range*eps_step) {
                    if (r0_j < 1.11 || r0_j > 2.10 || eps_j < 0.1 || eps_j > 5.0) {
                        continue;
                    }
                    if (r0_i == r0_j && eps_i == eps_j) {
                        continue;
                    }
                    
                    char r0stri[8], epsstri[8], r0strj[8], epsstrj[8];
                    sprintf(r0stri, "%1.2f", r0_i);
                    sprintf(epsstri, "%1.1f", eps_i);
                    sprintf(r0strj, "%1.2f", r0_j);
                    sprintf(epsstrj, "%1.1f", eps_j);
                    
                    total += matchfundist(descriptor[r0stri][epsstri], descriptor[r0strj][epsstrj]);
                    count ++;
                }
            }
            if (count > 0) {
                matrix << total/(double)count << "\t";
            }
            else {
                matrix << "0.0\t";            
            }
        }
        matrix << "\n";
    }
    matrix.close();
}