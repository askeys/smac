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
 *  Cell.h
 *  glotzilla
 *
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GLOTZILLA_STDCELL_H
#define _GLOTZILLA_STDCELL_H

#include <vector>
#include <map>

namespace smac
{
	/**
	@brief Contains a list of integers 
	@author Aaron Keys
	@ingroup opt
	
	Cell represents one cell in class CellList, which contains and manages many Cells
	*/
	class Cell
	{
		public:
			Cell();																///<Constructor
			virtual ~Cell();													///<Destructor
			
			virtual void clear();												///<Clears the cell memory space (contents, neighbors)
			virtual void clearContents();										///<Clears the contents of the cell, but not the neighbors

			virtual void insert(int);											///<Inserts an object at the back of the cell
			virtual void remove(int);											///<Removes an object from the cell

			void addNeighbor(int);												///<Add a neighboring cell to the list of neighbors

			int getNumberOfObjects();											///<Returns the number of objects in the cell
			int getNumberOfNeighbors() const;									///<Returns the number of neighboring cells
			std::vector<int> getNeighbors();									///<Returns a list of neighbors as stl vector<int>
			std::vector<int> getContents();									///<Returns a list of contents as stl vector<int>
			
		protected:
			std::vector<int> _neighbor;											///<List containing the neighboring cells of a cell
			std::vector<int> _tag;												///<A container to hold the tags of objects in the cell
			std::map<int, int> _index_of;										///Index of a given tag within the _tag array
	};
}
#endif

