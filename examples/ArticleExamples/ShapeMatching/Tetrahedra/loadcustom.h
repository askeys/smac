/*
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
*/#ifndef LOADDATAFILE
#define LOADDATAFILE

#include <smgz/smgz.h>
#include <fstream>
#include <cmath>

void loadcustom(const char* filename, int npart,
	smgz::coordlist_t& x, smgz::coordlist_t& x_cm, smgz::box_t& box)
{ 
	std::ifstream in_file (filename);
	assert(in_file.good());
	x.resize(npart*4);
	for (unsigned int i=0; i<x.size(); i++) {
		x[i] = smgz::coord_t(3, 0.0);
	}
	x_cm.resize(npart);
	for (unsigned int i=0; i<x_cm.size(); i++) {
		x_cm[i] = smgz::coord_t(3, 0.0);
	}
	
	double temp, q[4], R[3][3];
	double L;
	in_file >> L >> temp >> temp >> temp >> temp;
	
	box.period = std::vector<double> (3, L);
	box.periodic = std::vector<bool> (3, true);
	box.boxlo = std::vector<double> (3, -L/2.0);
	box.boxhi = std::vector<double> (3, L/2.0);
	
	double a_0[4][3];
	a_0 [0][0] = sqrt (2.0/3.0);
	a_0 [0][1] = 0;
	a_0 [0][2] = sqrt (1.0/3.0);
	a_0 [1][0] = -sqrt (2.0/3.0);
	a_0 [1][1] = 0;
	a_0 [1][2] = sqrt (1.0/3.0);
	a_0 [2][0] = 0;
	a_0 [2][1] = sqrt (2.0/3.0);
	a_0 [2][2] = -sqrt (1.0/3.0);
	a_0 [3][0] = 0;
	a_0 [3][1] = -sqrt (2.0/3.0);
	a_0 [3][2] = -sqrt (1.0/3.0);
	double sqrt_38 = sqrt(3.0/8.0);

	
	for (int i=0; i<npart; i++) {
		for (int j=0; j<3; j++) {
			in_file >> temp;
			x_cm[i][j] = temp;
		}
	
		for (int j=0; j<4; j++) {
			in_file >> q[j];
		}
 
		double q00 = q[0]*q[0];
		double q01 = q[0]*q[1];
		double q02 = q[0]*q[2];
		double q03 = q[0]*q[3];
		double q11 = q[1]*q[1];
		double q12 = q[1]*q[2];
		double q13 = q[1]*q[3];
		double q22 = q[2]*q[2];
		double q23 = q[2]*q[3];
		double q33 = q[3]*q[3];
 
		R[0][0] = q00 + q11 - q22 - q33;
		R[1][0] = 2.0 * ( q12 - q03 );
		R[2][0] = 2.0 * ( q02 + q13 );
		R[0][1] = 2.0 * ( q12 + q03 );
		R[1][1] = q00 - q11 + q22 - q33;
		R[2][1] = 2.0 * ( q23 - q01 );
		R[0][2] = 2.0 * ( q13 - q02 );
		R[1][2] = 2.0 * ( q01 + q23 );
		R[2][2] = q00 - q11 - q22 + q33;
		
		
		for (int j=0; j<4; j++) {
			for (int s=0; s<3; s++) {
				temp = 0;
				for (int t=0; t<3; t++) {
					temp +=  R[s][t] * a_0[j][t];
				}
				x[i*4+j][s] = temp*sqrt_38 + x_cm[i][s];
			}
		} 
	}
	
	in_file.close();
}
 
#endif
