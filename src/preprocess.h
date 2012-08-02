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
*//*
 *  f_coord.h
 *  glotzilla
 *
 *  Created by askeys on 9/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PREPROCESS
#define PREPROCESS

#include "typedefs.h"
#include "CellList.h"

namespace smac
{
	void perturb(coordlist_t& x, double max);
    void blur(coordlist_t& x, double max, int npts);
	
    //void radialsubset(coordlist_t& x, double rmin, double rmax );
	coordlist_t subset(coordlist_t&, clstr_t&);
    
	void unmap(coordlist_t& x, coord_t& x_c, box_t& box);
	void unmapcentroid(coordlist_t& x, box_t& box);

	coord_t centroid(coordlist_t& x, box_t& box);

	coord_t centroid(coordlist_t& x);
		
	void translate(coordlist_t& x, const coord_t& x_c);
	void untranslate(coordlist_t& x, const coord_t& x_c);

    struct extractsurface_info
    {
        extractsurface_info() : rcentroidcut(-1.0), range(-1.0), minnbrs(1), dim(3) {}
        double rcentroidcut;
        double range;
        int minnbrs;
        int dim;
    };
    void extractsurface(shapedata_t& s, box_t& box, std::vector<int>& surface, 
        extractsurface_info& info);
	
    /**
    \brief rotate a coordinate list by an angle (radians) about the x, y, z axes
    */
	void rotate(coordlist_t& x, double rotx, double roty, double rotz);

    /**
    \brief rotate a coordinate list by an angle (degrees) about the x, y, z axes
    */
	void rotated(coordlist_t& x, double rotx, double roty, double rotz);

    /**
    \brief rotate a coordinate list by an angle (radians) about a specified axis
    */
	void rotateaa(coordlist_t& x, double angle, coord_t& axis);

    /**
    \brief rotate a coordinate list by a specified quaternion
    */
	void rotateq(coordlist_t& x, std::vector<double> Q);
	
	void coordtovoxel(coordlist_t& x, CellList& voxels, int resolution=10);

	void transcentroid(coordlist_t& x);
	void transorigin(coordlist_t& x);
	
    /**
    \brief normalize a set of coordinates to the surface of a unit sphere
    */
    void normusphere(coordlist_t& x);
    /**
    \brief normalize a set of coordinates to lie within a unit box
    */
	void normubox(coordlist_t& x, double);

    /**
    \brief normalize a set of coordinates to lie within a unit box
    */
	void normubox(coordlist_t& x);

    /**
    \brief normalize a set of coordinates in a box of a specified size 
        to lie within a unit box
    */
	void normuball(coordlist_t& x);

    /**
    \brief generate a regual polygon of a specified size in the xy plane
    */
	void polygon(int, coordlist_t& x);

    /**
    \brief rescale coorindates by a scaling factor
    */
	void rescale(coordlist_t&x, double s);

    /**
    \brief convert 2d data to 3d data by adding zeros to the z dimension
    */
	void pad3d(coordlist_t& x);
	
    /**
    \brief checks whether each coordinate is within the specified box
    */
	void insidebox(coordlist_t&, box_t&, std::vector<bool>&);

    /**
    \brief checks whether each coordinate is within the specified radius from 
        the origin
    */
	void insidesphere(coordlist_t&, coord_t&, double, std::vector<bool>&);
	
}

#endif
