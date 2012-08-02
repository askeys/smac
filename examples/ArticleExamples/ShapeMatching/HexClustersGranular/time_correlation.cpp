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
#include <vmdstream/vmdstream.h>

using namespace std;
using namespace smac;

#define R0 1.65
#define R1 1.65

#define BONDCUT 0.99
#define FBONDCUT 0.5

void print_average(map<int, pair<double, int> >& totalcount, ofstream& f)
{
    for (map<int, pair<double, int> >::iterator i=totalcount.begin();
        i!=totalcount.end(); ++i) {
        int timewindow = i->first;
        pair<double, int> p = i->second;
        double total = p.first;
        double count = p.second;
        f << timewindow << "\t" << total/count << "\n";
    }
}


int main(int argc, char** argv)
{		
    string datapath = string(_SMAC_SRC_PATH) + "/example_data/2dliquid/";

    double density = 0.75;
    double box_size = sqrt(1024/density)*0.5;
	box_t box;
	box.period = vector<double>(3, 2*box_size);
	box.boxhi = vector<double>(3, box_size);
	box.boxlo = vector<double>(3, -box_size);
	box.boxhi[2] = 3.0;
	box.boxlo[2] = -3.0;
	box.period[2] = 6;
	box.periodic[2] = false;

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
	
    //set up some parameters for the localorder function
    localorder_info lo_info;
    lo_info.origin = CENTER_PARTICLE;

    //get a list of datafiles:
    vector<string> file_list;
    ostringstream ls_command;
    ls_command << "ls " << datapath << " | sort -n > files.txt";
    system(ls_command.str().c_str());
    ifstream f("files.txt");
    
    while (1) {
        string tmp;
        f >> tmp;
        if (f.fail()) {
            break;
        }
        file_list.push_back(tmp);
    }
    f.close();
    
    map<int, pair<double, int> > dot_liq;
    map<int, pair<double, int> > dot_hex;
    
    for (unsigned int i=0; i<file_list.size(); i+=100) {
        string file = datapath + file_list[i];
        xyzfile xyz;
        loadxyz(file.c_str(), xyz);
        coordlist_t x0 = xyz.x;
        for (int ii=0; ii<x0.size(); ii++) {
            pbc(x0[ii], box.period, box.periodic);
        }
        //find the crystal grains:
        clstrlist_t crystal_grains;
        crystalgrains(x0, box, c_info, crystal_grains);
        clstrlist_t clusters;
        shpdesclist_t sd_all;
        localorder(x0, box, R0, rmsdesc, &sd_args, clusters, sd_all, lo_info);
        vector<bool> in_crystal(1024, false);
        for (unsigned int ii=0; ii<crystal_grains.size(); ii++) {
            if (crystal_grains[ii].size() > 3) {
                for (unsigned int jj=0; jj<crystal_grains[ii].size(); jj++) {
                    in_crystal[crystal_grains[ii][jj]] = true;
                }
            }
        }
        
        for (unsigned int j=i; j<file_list.size(); j++) {
            cerr << file_list[i] << "\t" << file_list[j] << "\n"; 
            string file = datapath + file_list[j];
            xyzfile xyz;
            loadxyz(file.c_str(), xyz);
            coordlist_t xt = xyz.x;
            for (int ii=0; ii<xt.size(); ii++) {
                pbc(xt[ii], box.period, box.periodic);
            }
            clstrlist_t clusters;
            shpdesclist_t sd_all_t;
            localorder(xt, box, R0, rmsdesc, &sd_args, clusters, sd_all_t, lo_info);
            
            int delta = j-i;
            for (int ii=0; ii<sd_all_t.size(); ii++) {
                rmsfastassign(sd_all[ii], sd_all_t[ii]);
                if (sd_all[ii].size() > sd_all_t[ii].size()) {
                    sd_all[ii].resize(sd_all_t[ii].size());
                }
                if (sd_all_t[ii].size() > sd_all[ii].size()) {
                    sd_all[ii].resize(sd_all_t[ii].size());
                }
                
                double d = matchfundot(sd_all[ii], sd_all_t[ii]);
                if (in_crystal[ii]) {
                    if (dot_hex.find(delta) == dot_hex.end()) {
                        dot_hex[delta] = pair<double, int>(0.0, 0);
                    }
                    dot_hex[delta].first += d;
                    dot_hex[delta].second ++;
                }
                else {
                    if (dot_liq.find(delta) == dot_liq.end()) {
                        dot_liq[delta] = pair<double, int>(0.0, 0);
                    }
                    dot_liq[delta].first += d;
                    dot_liq[delta].second ++;                
                }
            }
        }
    }
    ofstream os1("dot_liq.txt");
    ofstream os2("dot_hex.txt");
    print_average(dot_liq, os1);
    print_average(dot_hex, os2);
    os1.close();
    os2.close();
}