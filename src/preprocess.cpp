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
 *  f_coord.cpp
 *  glotzilla
 *
 *  Created by askeys on 9/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "preprocess.h"
#include "space.h"
#include "rand.h"
#include "cluster.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace smac
{	
	void perturb(coordlist_t& x, double max)
	{
		for (unsigned int i=0; i<x.size(); i++) {
			x[i][0] += randgauss()*max;
			x[i][1] += randgauss()*max;
			x[i][2] += randgauss()*max;		
		}
	}
	
	void blur(coordlist_t& x, double max, int npts)
	{
		if (x.size() == 0) {
			return;
		}
		
		coordlist_t x_rot;
        //x_rot.reserve(x.size());
		int dim = x[0].size();		
		for (unsigned int i=0; i<x.size(); i++) {
			for (int j=0; j<npts; j++) {
				coord_t xnew = x[i];
				for (int k=0; k<dim; k++) {
					xnew[k] += randgauss()*max;
				}
				x_rot.push_back(xnew);
			}
		}
		x.insert(x.end(), x_rot.begin(), x_rot.end());
	}

	void subset_radial(coordlist_t& x, 
		double rmin, 
		double rmax )
	{
		double rminsq = rmin*rmin;
		double rmaxsq = rmax*rmax;
		coordlist_t x_subset;
		coord_t o(3, 0.0);
	
		for (unsigned int i=0; i<x.size(); i++) {
			double rsq = distancesq(x[i], o);
			if (rsq >= rminsq && rsq < rmaxsq) {
				x_subset.push_back(x[i]);
			}
		}
		x = x_subset;
	}
    
    /**
    \param x is the coordinates to unmap
    \param box is the simulation box
    
    Unmaps the periodic boundary conditions and places the object
    at the centroid
    */
	void unmapcentroid( coordlist_t& x, box_t& box)
	{
		coord_t x_c = centroid(x, box);
		unmap(x, x_c, box);
		x_c[0] = -x_c[0]; x_c[1] = -x_c[1]; x_c[2] = -x_c[2];
		translate(x, x_c);
	}
	
  /**
    \param x is the coordinates to unmap
    \param box is the simulation box
    
    Unmaps the periodic boundary conditions and places the object
    at the a specified coordinate
    */
	void unmap(coordlist_t& x, coord_t& x_c, box_t& box)
	{		
		if (x.size() == 0) {
			return;
		}
		int dim = x[0].size();
		
		for (unsigned int i=0; i<x.size(); i++) {
			for (int k=0; k<dim; k++) {
				double dx = x[i][k] - x_c[k];
				pbc(dx, box.period[k], box.periodic[k]);
				x[i][k] = dx;
			}
		}
	}
	
	coord_t centroid(coordlist_t& x)
	{
		int n = x.size();
		if (n == 0) {
			return coord_t(3, 0.0);
		}
		int dim = x[0].size();
		
		coord_t xc(dim, 0.0);
		
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				xc[j] += x[i][j];
			}
		}

		for (int j=0; j<dim; j++) {
			xc[j] /= (double)n;
		}
		return xc;		
	}

	coord_t centroid(coordlist_t& x, box_t& box)
	{
		assert(x.size() > 0);
		coord_t rc(3, 0.0);
		coord_t x_ref = x[0];
		
		for (unsigned int i=0; i<x.size(); i++) {
			coord_t xi(3);
			xi[0] = x[i][0] - x_ref[0];
			xi[1] = x[i][1] - x_ref[1];
			xi[2] = x[i][2] - x_ref[2];
			pbc(xi, box.period, box.periodic);
			rc[0] += xi[0]; 
			rc[1] += xi[1]; 
			rc[2] += xi[2];
		}
		
		rc[0] /= (double)x.size(); 
		rc[1] /= (double)x.size(); 
		rc[2] /= (double)x.size();
		
		rc[0] += x_ref[0];
		rc[1] += x_ref[1];
		rc[2] += x_ref[2];
		
		return rc;
	}
	
	void translate(coordlist_t& x, const coord_t& x_c)
	{
		for (unsigned int i=0; i<x.size(); i++) {
			x[i][0] += x_c[0]; x[i][1] += x_c[1]; x[i][2] += x_c[2];
		}
	}

	void untranslate(coordlist_t& x, const coord_t& x_c)
	{
		for (unsigned int i=0; i<x.size(); i++) {
			x[i][0] -= x_c[0]; x[i][1] -= x_c[1]; x[i][2] -= x_c[2];
		}
	}
	
    /**
    \param x is the list of coordinates to normalize
    */
	void normusphere(coordlist_t& x)
	{
		for (unsigned int i=0; i<x.size(); i++) {
			double r = sqrt(x[i][0]*x[i][0]+x[i][1]*x[i][1]+x[i][2]*x[i][2]);
			if (r!=0.0) {
				x[i][0]/=r; x[i][1]/=r; x[i][2]/=r;
			}
		}
	}
	
	void rotate(coordlist_t& x,
		double rotx, double roty, double rotz)
	{
		if (rotx != 0.0) {
			coord_t xaxis(3, 0.0);
			xaxis[0] = 1.0;
			rotateaa(x, rotx, xaxis);
		}

		if (roty != 0.0) {
			coord_t yaxis(3, 0.0);
			yaxis[0] = 1.0;
			rotateaa(x, roty, yaxis);
		}
		
		if (rotz != 0.0) {
			coord_t zaxis(3, 0.0);
			zaxis[2] = 1.0;
			rotateaa(x, rotz, zaxis);
		}
	}

	void rotated(coordlist_t& x,
		double rotx, double roty, double rotz)
	{
		rotate(x, 
			rotx * M_PI / 180.0, roty * M_PI / 180.0, rotz * M_PI / 180.0);
	}

	void rotateaa(coordlist_t& x,
		double angle, coord_t& axis)
	{
		assert(axis.size() >= 3);
		double norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
		axis[0] /= norm; axis[1] /= norm; axis[2] /= norm;

		coord_t q(4, 0.0);
		q[0] = cos(angle/2.0);
		q[1] = axis[0] * sin(angle/2.0);
		q[2] = axis[1] * sin(angle/2.0);
		q[3] = axis[2] * sin(angle/2.0);
		rotateq(x, q);
	}
	
	void rotateq(coordlist_t& x, std::vector<double> Q)
	{	
		double norm = sqrt(Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]);
		Q[0] /= norm; Q[1] /= norm; Q[2] /= norm; Q[3] /= norm;
		
		double AXX = Q[0] * Q[0] + Q[1] * Q[1] - Q[2] * Q[2] - Q[3] * Q[3];
		double AXY = 2.0 * ( Q[1] * Q[2] + Q[0] * Q[3] );
		double AXZ = 2.0 * ( Q[1] * Q[3] - Q[0] * Q[2] );
		double AYX = 2.0 * ( Q[1] * Q[2] - Q[0] * Q[3] );
		double AYY = Q[0] * Q[0] - Q[1] * Q[1] + Q[2] * Q[2] - Q[3] * Q[3];
		double AYZ = 2.0 * ( Q[2] * Q[3] + Q[0] * Q[1] );
		double AZX = 2.0 * ( Q[1] * Q[3] + Q[0] * Q[2] );
		double AZY = 2.0 * ( Q[2] * Q[3] - Q[0] * Q[1] );
		double AZZ = Q[0] * Q[0] - Q[1] * Q[1] - Q[2] * Q[2] + Q[3] * Q[3];
		
		int n = x.size();
		for (int i=0; i<n; i++) {
			//create vector between a vertex and the center
			
			norm = sqrt(
				x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2]);
			
			double unit_vector[3];
			if (norm == 0.0) {
				unit_vector[0] = unit_vector[i] = unit_vector[2] = 0.0;
			}
			else {
				unit_vector[0] = x[i][0] / norm; 
				unit_vector[1] = x[i][1] / norm; 
				unit_vector[2] = x[i][2] / norm;
			}
				
			double temp[3];
			temp[0] = AXX * unit_vector[0] + 
				AXY * unit_vector[1] + AXZ * unit_vector[2];
			temp[1] = AYX * unit_vector[0] + 
				AYY * unit_vector[1] + AYZ * unit_vector[2];
			temp[2] = AZX * unit_vector[0] + 
				AZY * unit_vector[1] + AZZ * unit_vector[2];
		
			x[i][0] = temp[0] * norm;
			x[i][1] = temp[1] * norm;
			x[i][2] = temp[2] * norm;
			
		}
	}
	
	void coordtovoxel(coordlist_t& x, CellList& cell_list, int resolution)
	{
		cell_list.clear();
		if (x.size() == 0) {
			return;
		}
		
		if (x[0].size() != 3) {
			std::cerr << "coordtovoxel: Error. Data must be 3 dimensional\n";
			return;
		}
		
		coord_t xmin(3, 1e10), xmax(3, -1e10);
		for (unsigned int i=0; i<x.size(); i++) {
			for (int j=0; j<3; j++) {
				xmin[j] = std::min(xmin[j], x[i][j]);
				xmax[j] = std::max(xmax[j], x[i][j]);
			}
		}
		
		coord_t box(3);
		for (int j=0; j<3; j++) {
			box[j] = xmax[j] - xmin[j];
		}
		double lmax = *(std::max_element(box.begin(), box.end()));
		for (int j=0; j<3; j++) {
			box[j] = 1.1*lmax;
		}
		coord_t boxlo(3), boxhi(3);
		for (int j=0; j<3; j++) {
			boxlo[j] = 0.5*(xmin[j] + xmax[j]) - 0.5*box[j];
			boxhi[j] = 0.5*(xmin[j] + xmax[j]) + 0.5*box[j];
		}
		
		cell_list.setInteractionRange(lmax/(double)resolution);
		cell_list.setBox(boxlo, boxhi);
		
		for (unsigned int i=0; i<x.size(); i++) {
			cell_list.insert(i, x[i]);
		}
	}
	
	void transcentroid(coordlist_t& x)
	{
		int n = x.size();
		if (n == 0) {
			return;
		}
		int dim = x[0].size();
		
		coord_t cm = centroid(x);

		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				x[i][j] -= cm[j];
			}
		}
	}

	void transorigin(coordlist_t& x)
	{
		int n = x.size();
		if (n == 0) {
			return;
		}
		int dim = x[0].size();

		coord_t xmin(3, 1e32);
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				xmin[j] = std::min(xmin[j], x[i][j]);
			}
		}
		
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				x[i][j] -= xmin[j];
			}
		}
		
	}

    /**
    \param x is the list of coordinates to normalize
    \param boxdim is the size of the original box
    */
	void normubox(coordlist_t& x, double boxdim)
	{
		int n = x.size();
		if (n == 0) {
			return;
		}
		transorigin(x);
		int dim = x[0].size();
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				x[i][j] /= boxdim;
				x[i][j] -= 0.5;
				//assert(fabs(x[i][j]) < 0.5);
			}
		}
	}

    /**
    \param x is the list of coordinates to normalize
    */
	void normubox(coordlist_t& x)
	{
		int n = x.size();
		if (n == 0) {
			return;
		}
		
		transorigin(x);

		int dim = x[0].size();
		double xmax = -1e30;
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				xmax = std::max(xmax, x[i][j]);
			}
		}
		double boxdim = xmax*1.01;
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				x[i][j] /= boxdim;
				x[i][j] -= 0.5;
				//assert(fabs(x[i][j]) < 0.5);
			}
		}
		//transcentroid(x);
	}
	
    void extractsurface(
        shapedata_t& s, box_t& box, std::vector<int>& surface, 
        extractsurface_info& info)
    {
        if (info.range <= 0.0) {
            std::cerr << "ERROR: extractsurface: argument info.range";
            std::cerr << "must be specified and set >= 0.0\n";
            throw std::runtime_error("ERROR setting extractsurface_info.range");
        }
        if (info.rcentroidcut <= 0.0) {
            std::cerr << "ERROR: extractsurface: argument info.rcentroidcut";
            std::cerr << "must be specified and set >= 0.0\n";
            throw std::runtime_error("ERROR setting extractsurface_info.rcentroidcut");
        }

        //1) get neighbors
        //2) count neighbors
        //3) get centroid
        // a particle on the surface will have ~0.5x as many neighbors as max
        // a particle on the surface will have lobsided cluster
        
		clstrrulerange_info args;
		args.x = s.x;
		args.box = box;
		args.dim = info.dim;
		args.rcutsq = info.range*info.range;
		
		clstrlist_t cluster;
		clstrshortrange(s.x, box, info.range, clstrrulerange, &args, cluster);
        
        double centroidcutsq = info.rcentroidcut*info.rcentroidcut;
        
        for (unsigned int i=0; i<cluster.size(); i++) {
            if (cluster[i].size() < info.minnbrs) {
                continue;
            }
			coordlist_t x_cluster(cluster[i].size());
			for (unsigned int j=0; j<cluster[i].size(); j++) {
				x_cluster[j] = s.x[cluster[i][j]];
			}
            coord_t x_c = s.x[cluster[i][0]];
            unmap(x_cluster, x_c, box);
            x_c = centroid(x_cluster);
            double rsq = x_c[0]*x_c[0] + x_c[1]*x_c[1] + x_c[2]*x_c[2];
            if (rsq > centroidcutsq) {
                surface.push_back(i);
            }
        }
    }

	/**
    \param x is the coordinates to normalize
	*/
	void normuball(coordlist_t& x)
	{
        double rsqmax = 0.0;
        for (unsigned int i=0; i<x.size(); i++) {
			double rsq = x[i][0]*x[i][0]+x[i][1]*x[i][1]+x[i][2]*x[i][2];
            rsqmax = std::max(rsq, rsqmax);
            //std::cerr << sqrt(rsq) << std::endl;
        }
        if (rsqmax <= 0.0) {
            return;
        }
        double invrmax = 1.0/sqrt(rsqmax);
        //std::cerr << "max: " << sqrt(rsqmax) << std::endl;
		for (unsigned int i=0; i<x.size(); i++) {
            x[i][0]*=invrmax; x[i][1]*=invrmax; x[i][2]*=invrmax;
		}

        /*
		int n = x.size();
		if (n == 0) {
			return;
		}
		int dim = x[0].size();
		double rmax = -1;
		for (int i=0; i<n; i++) {
			double r = 0.0;
			for (int j=0; j<dim; j++) {
				r += x[i][j]*x[i][j];
			}
			rmax = std::max(r, rmax);
		}
        double invrmax = 1.0/sqrt(rmax);
		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				x[i][j] *= invrmax;
			}
		}
        */
	}
	
	void polygon(int n, coordlist_t& x)
	{		
		x.clear();
		double step = 2*M_PI / (double)n;
		for (int i=0; i<n; i++) {
			coord_t xi(2);
			xi[0] = cos(i*step);
			xi[1] = sin(i*step);
			x.push_back(xi);
		}
	}
	
	void rescale(coordlist_t&x, double s)
	{
		int n = x.size();
		if (n == 0) {
			return;
		}
		int dim = x[0].size();

		for (int i=0; i<n; i++) {
			for (int j=0; j<dim; j++) {
				x[i][j] *= s;
			}
		}
	}
	
	void pad3d(coordlist_t& x)
	{
		int n = x.size();
		if (n == 0) {
			return;
		}
		int dim = x[0].size();
		
		coordlist_t xt = x;
		x.clear();
		for (int i=0; i<n; i++) {
			coord_t xi(3, 0.0);
			for (int j=0; j<dim; j++) {
				xi[j] = xt[i][j];
			}
			x.push_back(xi);
		}
	}
	
	void voxelize(coordlist_t&)
	{
		
		/*
		if (boxsize <= 0.0) {
			normubox(x);
		}
		else {
			normubox(x, boxsize);
		}

		int n = x.size();
		int ncell = pow(resolution, dim);
		
		coord_t factor(3, 1.0);
		if (dim > 1) {
			factor[1] = resolution;
		}
		if (dim > 2) {
			factor[2] = resolution*resolution;
		}
		
		std::vector<int> val(ncell, 0.0);
		for (int i=0; i<n; i++) {
			int cell_id = 0;
			for (int j=0; j<dim; j++) {
				cell_id += factor[j]*((int)((x[i][j]+0.5)*resolution));
			}
			//std::cerr << x[i][0] <<"\t" << x[i][1] <<"\t"<< cell_id <<"\n";
			val[cell_id] ++;
		}
		*/
	}
	
	void insidebox(coordlist_t& x, box_t& box, std::vector<bool>& inside)
	{
		inside.resize(x.size(), true);
		for (unsigned int i=0; i<x.size(); i++) {
			for (int k=0; k<3; k++) {
				if (x[i][k] < box.boxlo[k] || x[i][k] >= box.boxhi[k]) {
					inside[i] = false;
					break;
				}
			}
		}
	}
		
	void insidesphere(
		coordlist_t& x, coord_t& xc, double r, std::vector<bool>& inside)
	{
		double cutsq = r*r;
		inside.resize(x.size(), true);
		for (unsigned int i=0; i<x.size(); i++) {
			if (distancesq(x[i], xc) > cutsq) {			
				inside[i] = false;
			}			
		}
	}
    
    coordlist_t subset(coordlist_t& x, clstr_t& index)
    {
        coordlist_t xs(index.size());
        for (unsigned int i=0; i<index.size(); i++) {
            xs[i] = x[index[i]];
        }
        return xs;
    }

}
