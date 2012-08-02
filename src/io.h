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
 *  io.h
 *  libsmac
 *
 *  Created by askeys on 12/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SMIO
#define SMIO

#include "typedefs.h"
#include <iostream>
#include <string>

namespace smac 
{
	void load(const char*, coordlist_t&, arg_t=0x0);
	void loadxyz(const char* filename, coordlist_t& x, arg_t=0x0);
	void loadlmptraj(const char* filename, coordlist_t& x, arg_t);
	void loadlmptraj( std::istream& is, coordlist_t& x, arg_t arg);
	struct xyz_info 
	{
		xyz_info() : outstream(0x0), instream(0x0) 
		{ 
			xres = coord_t(3, 0.0);
		}
		
		box_t box;
		std::vector<std::string> type;
		std::ostream* outstream;
		std::istream* instream;
		std::map< std::string, int > resevoir;
		coord_t xres;
	};
	
	struct lmptraj_info
	{
		lmptraj_info() {}		
		box_t box;
		std::vector<std::string> type;
		std::istream* instream;
	};

    /**
    \brief datastructure to hold rasmol xyz data
    \ingroup io
    */
    struct xyzfile
    {
        unsigned int n;
        coordlist_t x;
        std::string commentstr;
        std::vector<std::string> type;
    };
    
    void loadxyz(const char*, xyzfile&);
    void loadxyz(std::istream&, xyzfile&);

    /**
    \brief datastructure to hold lammpstrj data
    \ingroup io
    */
    struct lammpstrjfile
    {
        unsigned int n;
        coordlist_t x;
        std::vector<int> type;
        int timestep;
        box_t box;
    };
    
    void loadlammpstrj(const char*, lammpstrjfile&);
    void loadlammpstrj(std::istream&, lammpstrjfile&);


	void savexyz(const char* filename, coordlist_t&x, xyz_info& info);
	void savevarxyz(const char* filename, coordlist_t&, xyz_info&);

	std::ostream& operator << (std::ostream&, box_t&);	
	std::ostream& operator << (std::ostream&, function1d_t&);	
	std::ostream& operator << (std::ostream&, coordlist_t&);
	std::ostream& operator << (std::ostream&, coord_t&);
	std::ostream& operator << (std::ostream&, shpdesc_t&);

	std::istream& operator >> (std::istream&, coordlist_t&);
	std::istream& operator >> (std::istream&, coord_t&);
	std::istream& operator >> (std::istream&, shpdesc_t&);

	void deltype(coordlist_t& x, std::vector<std::string>& types, std::string);
}

#endif
