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
 *  voxelize.cpp
 *  libsmac
 *
 *  Created by askeys on 2/8/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "voxel.h"
#include "preprocess.h"
#include <iostream>

#ifdef _USE_IPGZ
#include "../3rdparty/ipgz/IPGZ.h"
#endif


#include <assert.h>
#include <cstdlib>

namespace smac
{
	void voxelize(coordlist_t& x, voxel_t& voxel, int nx, int ny, int nz,
		coord_t box)
	{
		int n = x.size();
		assert(n > 0);

		coordlist_t xcopy = x;
		if(box.size() ==0) {
			normubox(xcopy);
		}
		else {
			double boxl = -1;
			for (unsigned int i=0; i<box.size(); i++) {
				boxl = std::max(box[i], boxl)*1.01;
			}
			normubox(xcopy, boxl);
		}

		int dim = x[0].size();
		int factor[3], ncell[3];
		factor[0] = ny*nz; factor[1] = nz; factor[2] = 1;
		ncell[0] = nx; ncell[1] = ny; ncell[2] = nz;

		voxel.resize(nx*ny*nz);
		for (unsigned int i=0; i<voxel.size(); i++) {
			voxel[i] = 0.0;
		}

		for (int i=0; i<n; i++) {
			int index = 0;
			for (int j=0; j<dim; j++) {
				index += factor[j]*(int)((xcopy[i][j]+0.5)*ncell[j]);
			}
			voxel[index] ++;
		}
	}

	void voxnorm(voxel_t& voxel)
	{
		double sum = 0.0;
		for (unsigned int i=0; i<voxel.size(); i++) {
			sum += voxel[i];
		}
		for (unsigned int i=0; i<voxel.size(); i++) {
			voxel[i] /= sum;
		}
	}

	void voxfiltergauss3d(voxel_t& voxel, int nx, int ny, int nz, int f)
	{
#ifdef _USE_IPGZ
		if (f==-1) {
			f = nx/8;
		}

		int nzny = nz*ny;
		IPGZ::PixelArray3D<float> pixels(nx, ny, nz), pixels_filtered;

		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				for (int k=0; k<nz; k++) {
					int index = i*nzny + j*ny + k;
					pixels.setValue(i, j, k, voxel[index]);
				}
			}
		}
		IPGZ:: gaussianFilter3DfastPBC(&pixels, &pixels_filtered, f, f, f);
		for (int i=0; i<nx; i++) {
			for (int j=0; j<ny; j++) {
				for (int k=0; k<nz; k++) {
					int index = i*nzny + j*ny + k;
					voxel[index] = pixels_filtered.getValue(i, j, k);
				}
			}
		}
#endif
	}

	void voxfilterfft3d(voxel_t& voxel, int nx, int ny, int nz)
	{
		int nzny = nz*ny;
		double twopi = 2.0*M_PI;
		double halfcellx = (2.0*M_PI)/(double)nx/2.0;
		double halfcelly = (2.0*M_PI)/(double)ny/2.0;
		double halfcellz = (2.0*M_PI)/(double)nz/2.0;
		std::vector< std::complex<double> > ft(nx*ny*nz);

		for (unsigned int i=0; i<ft.size(); i++) {
			ft[i] = std::complex<double>(0.0,0.0);
		}

		for (int q0=0; q0<nx; q0++) {
		for (int q1=0; q1<ny; q1++) {
		for (int q2=0; q2<nz; q2++) {

		for (int i=0; i<nx; i++) {
			double x = i/(double)nx*twopi + halfcellx;
			for (int j=0; j<ny; j++) {
				double y = j/(double)ny*twopi + halfcelly;
				for (int k=0; k<nz; k++) {
					double z = k/(double)nz*twopi + halfcellz;
					int index = i*nzny + j*ny + k;
					double real = voxel[index]*cos(x*q0 + y*q1 + z*q2);
					double imag = voxel[index]*sin(x*q0 + y*q1 + z*q2);
					ft[index] += std::complex<double>(real, imag);
				}
			}
		}
		}
		}
		}

		std::vector<double> ftabs(nx*ny*nz, 0);

		for (int q0=0; q0<nx; q0++) {
			for (int q1=0; q1<ny; q1++) {
				for (int q2=0; q2<nz; q2++) {
					int index = q0*nzny + q1*ny + q2;
					ftabs[index] = ft[index].real()*ft[index].real() +
									ft[index].imag()*ft[index].imag();
				}
			}
		}

		for (unsigned int i=0; i<voxel.size(); i++) {
			voxel[i] = 0.0;
		}

		for (int q0=0; q0<nx; q0++) {
		for (int q1=0; q1<ny; q1++) {
		for (int q2=0; q2<nz; q2++) {

		for (int i=0; i<nx; i++) {
			double x = i/(double)nx*twopi + halfcellx;
			for (int j=0; j<ny; j++) {
				double y = j/(double)ny*twopi + halfcelly;
				for (int k=0; k<nz; k++) {
					double z = k/(double)nz*twopi + halfcellz;
					int index = i*nzny + j*ny + k;
					voxel[index] += ftabs[index]*exp(x*q0 + y*q1 + z*q2);
				}
			}
		}
		}
		}
		}
		//voxnorm(voxel);
	}
}
