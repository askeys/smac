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
 *  pentagonal_dipyramids.h
 *  smgz
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

// The point of this example is to demonstrate matching when disordered 
// structres are possible.  This code is similar (although not identical) to 
// the one used to create the data for Fig 2, Nature.

#include <smgz/smgz.h>

#include "loadcustom.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>

using namespace std;
using namespace smgz;

//a naive function to build up a polygon cluster
//(can be made much faster by pruning paths that obviously will not form a 
//closed polygon)
void branch(int start, int i, int n, 
	clstrlist_t& cluster, clstrlist_t& polygon, 
	clstr_t& icluster, map<string, bool>& counted)
{	
	if (n!=1) {
		clstr_t iicluster = icluster;
		for (int j=1; j<cluster[i].size(); j++) {
			if (find(iicluster.begin(), iicluster.end(), cluster[i][j]) == iicluster.end()) {
				iicluster.push_back(cluster[i][j]);
				branch(start, cluster[i][j], n-1, cluster, polygon, iicluster, counted);
				iicluster.pop_back();
			}
		}	
	}
	else {
		for (int j=1; j<cluster[i].size(); j++) {
			if (cluster[i][j] == start) {
				clstr_t icluster_sorted = icluster;
				sort(icluster_sorted.begin(), icluster_sorted.end());
				ostringstream polystring;
				for (int k=0; k<icluster_sorted.size(); k++) {
					polystring << icluster_sorted[k] << ".";
				}				
				if(counted.find(polystring.str()) == counted.end()) {
					//std::cerr << polystring.str() << std::endl;
					counted[polystring.str()] = true;
					polygon.push_back(icluster);
				}
				return;
			}
		}
	}
}

//compute the edge directions for a pentagonal dipyramid
void edgedirections(coordlist_t& xtet, coordlist_t& xdir)
{
	xdir.clear();	
	coord_t xedge, xedgetemp(3);

	for (int i=0; i<5; i++) {
		double max = 0;
		for (int j=0; j<4; j++) {			
			for (int k=j+1; k<4; k++) {
				xedgetemp[0] = xtet[i*4 + k][0] + (xtet[i*4 + j][0] - xtet[i*4 + k][0]) / 2.0; 
				xedgetemp[1] = xtet[i*4 + k][1] + (xtet[i*4 + j][1] - xtet[i*4 + k][1]) / 2.0; 
				xedgetemp[2] = xtet[i*4 + k][2] + (xtet[i*4 + j][2] - xtet[i*4 + k][2]) / 2.0;
				//cerr << xedgetemp[0] <<"\t"<< xedgetemp[1] <<"\t"<< xedgetemp[2] <<"\n";
				double dist = xedgetemp[0]*xedgetemp[0] + xedgetemp[1]*xedgetemp[1] + xedgetemp[2]*xedgetemp[2];
				if ( dist > max) {
					max = dist;
					xedge = xedgetemp;
				}
			}
		}
		xdir.push_back(xedge);
	}
}

int main(int argc, char** argv)
{	
	//create the reference structure (just a regular pentagon):
	coordlist_t xref(5);
	double twopiover5 = 2.0*M_PI/5.0;
	for (int i=0; i<5; i++) {
		coord_t xrefi(3, 0.0);
		xrefi[0] = cos(i*twopiover5);
		xrefi[1] = sin(i*twopiover5);
		xref[i] = xrefi;
	}
	
	//create the reference shape descriptor:
	shpdesc_t ref_descriptor;
	fourier_info fargs;
	fargs.invariant = INVARIANT_Q;
	fargs.frequency.clear();
	for (int i=0; i<=5; i++) {
		fargs.frequency.push_back(i+5);
	}	
    
    shapedata_t shaperef(xref);
	fourierdesc3d(shaperef, ref_descriptor, &fargs);
	
	//get the input data:	
	coordlist_t x, xcm;
	box_t box;	

	int npart = atoi(argv[1]);
	loadcustom(argv[2], npart, x, xcm, box);
	double range = atof(argv[3]);
	
	clstrrulerange_info cargs;
	cargs.x = xcm;
	cargs.box = box;
	cargs.rcutsq = range*range;
	clstrlist_t cluster;
	double matchcut = atof(argv[4]);
	
	vector<int> npoly(npart, 0);
	vector<bool> inpoly(npart, false);
	int totalpoly = 0;
	
	clstrshortrange(xcm, box, range, clstrrulerange, &cargs, cluster);
	
	map<string, bool> counted;
	clstrlist_t polygon;
	for (int i=0; i<npart; i++) {
		clstr_t icluster(1, i);
		branch(i, i, 5, cluster, polygon, icluster, counted);
	}	
	
	for (int i=0; i<polygon.size(); i++) {
		coordlist_t xtemp;
		for (int j=0; j<polygon[i].size(); j++) {
			int idx = polygon[i][j];
			for (int k=0; k<4; k++) {
				coord_t xidx(3);
				xidx[0] = x[idx*4 + k][0];
				xidx[1] = x[idx*4 + k][1];
				xidx[2] = x[idx*4 + k][2];
				xtemp.push_back(xidx);
			}
		}
		
		coordlist_t xpolygon;
		
		unmapcentroid(xtemp, box);
		transcentroid(xtemp);
		edgedirections(xtemp, xpolygon); 
		
		shpdesc_t descriptor;											
		transcentroid(xpolygon);

        shapedata_t shapepolygon(xpolygon);
		fourierdesc3d(shapepolygon, descriptor, &fargs);			
		double match = matchfundist(descriptor, ref_descriptor, 0, 1); 

#ifdef DEBUG
		cerr << "\nmatch: " << match << "\n\n";

		for (int j=0; j<xtemp.size(); j++) {
			cout << xtemp[j][0] <<"\t" << xtemp[j][1] <<"\t" << xtemp[j][2] <<"\n";
		}
		cerr << "\n";
		
		for (int k=0; k<xpolygon.size(); k++) {
			cerr << xpolygon[k][0] <<"\t"<< xpolygon[k][1] <<"\t"<< xpolygon[k][2] <<"\n";
		}
#endif	
		
		if (match > matchcut) {
			totalpoly++;
			for (int k=0; k<polygon[i].size(); k++) {
				inpoly[polygon[i][k]] = true;
				npoly[polygon[i][k]] ++;
			}
		}
	}
	
	int ninpoly = 0, npolyper = 0;
	for (int i=0; i<npart; i++) {
		if (inpoly[i]) {
			ninpoly++;
		}
		npolyper += npoly[i];
	}
	cout << totalpoly << "\t";
	cout << npolyper / double(npart) << "\t";
	cout << ninpoly / double(npart) << "\n";
		
	return 0;
}