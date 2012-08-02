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
 *  typedefs.h
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SMTYPEDEFS
#define SMTYPEDEFS

#include <vector>
#include <complex>
#include <sstream>
#include <map>

#include "config.h"

namespace smac
{

	/**
	\typedef std::complex<double> component_t;
	\briefA component in a shape descriptor vector
	\ingroup types
	*/

	/**
	\typedef std::vector<component_t> shpdesc_t;
	\briefA shape descriptor vector
	\ingroup types
	*/
	/**
	\typedef std::vector<shpdesc_t> shpdesclist_t;
	\briefA list of shape descriptor vectors
	\ingroup types
	*/
	/**
	\typedef std::vector<double> coord_t;
	\briefA cartesian coordinate
	\ingroup types
	*/
	/**
	typedef std::vector<coord_t> coordlist_t;
	\briefA list of cartesian coordinates
	\ingroup types
	*/
	/**
	\typedef std::vector<bool> bool_array_t;
	\briefAn array of boolean values
	\ingroup types
	*/
	/**
	\typedef std::vector<double> voxel_t;
	\ingroup types
	*/
	/**
	\typedef std::vector<int> clstr_t;
	\ingroup types
	*/
	/**
	\typedef std::vector< clstr_t > clstrlist_t;
	\ingroup types
	*/
	/**
	\typedef std::map<double, double> function1d_t;

	\typedef void* arg_t;
	\ingroup types
	*/
	/**
	\typedef void (*shpdescfun_t) (shapedata_t&, shpdesc_t&, arg_t);
	\ingroup types
	*/
	/**
	\typedef void (*assignfun_t) (shpdesc_t&, shpdesc_t&, arg_t);
	\ingroup types
	*/
	/**
	\typedef double (*matchfun_t) (shpdesc_t&, shpdesc_t&, double, double);
	\ingroup types
	*/
	/**
	\typedef bool (*clstrrule_t)(int, int, arg_t);
	\ingroup types
	*/
	/**
	\typedef bool (*loadfun_t)(const char*, coordlist_t&, arg_t);
	\ingroup types
	\brief Pointer to a load function
	*/


	typedef std::complex<double> component_t;
	typedef std::vector<component_t> shpdesc_t;
	typedef std::vector<shpdesc_t> shpdesclist_t;
	typedef std::vector<double> coord_t;
	typedef std::vector< std::vector<double> > coordlist_t;
	typedef std::vector<bool> bool_array_t;
	typedef std::vector<double> voxel_t;
	typedef std::vector<int> clstr_t;
	typedef std::vector< clstr_t > clstrlist_t;
	typedef std::map<double, double> function1d_t;
	typedef std::vector<double> weight_t;

	/**
	\brief a struct containing shape data
	*/
	struct shapedata_t
	{
		shapedata_t()
            {}
		shapedata_t(coordlist_t& xi, weight_t& wi) :
			x(xi), f(wi) {}
		shapedata_t(coordlist_t& xi) :
			x(xi), f(weight_t(xi.size(), 1.0)) {}
		coordlist_t x;
		weight_t f;
	};

	typedef void* arg_t;
	typedef void (*shpdescfun_t) (shapedata_t&, shpdesc_t&, arg_t);
	typedef double (*matchfun_t) (shpdesc_t&, shpdesc_t&, double, double);
	typedef bool (*clstrrule_t)(int, int, arg_t);
	typedef void (*loadfun_t)(const char*, coordlist_t&, arg_t);
	typedef void (*assignfun_t)(shpdesc_t&, shpdesc_t&, arg_t);

	/**
	\brief a struct containing box information
	*/
	struct box_t
	{
		box_t() :
			boxlo(3,0),
			boxhi(3,0),
			period(3, 0),
			periodic(3, true) {}

		coord_t boxlo;
		coord_t boxhi;
		coord_t period;
		std::vector<bool> periodic;
	};

}

#endif
