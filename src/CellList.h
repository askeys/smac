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
 *  CellList.h
 *  glotzilla
 *
 *  Created by Aaron Keys on 9/30/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GLOTZILLA_STDCELLLIST_H
#define _GLOTZILLA_STDCELLLIST_H

#include "Cell.h"
#include <map>
#include <utility>
#include "typedefs.h"

namespace smac
{
	/**
	@brief Contains and manages many Cell(s).
	@author Aaron Keys
	@ingroup opt
	*/
	class CellList
	{
		public:
			CellList();															///<Constructor
			virtual ~CellList();												///<Destructor
			
			int getNumberOfCells() const;										///<Returns number of Cells in the CellList
			int getCellContaining(coord_t&) const;								///<Returns the cell containing a position in the box
			int getCellContaining(int);											///<Returns the cell containing an object

			coord_t getPositionOfCell(int);
			Cell getCell(int);													///<Returns cell at index i
			
			void setBox(coord_t&);												///<sets the box dimensions
			void setBox(coord_t&, coord_t&);									///<sets the box min, max
			void setBox(
				coord_t&, coord_t&, std::vector<bool>&);						///<sets the box min, max, periodicity
			void setMinimumCellSize(double);									///<sets the smallest possible cell dimension
			void setInteractionRange(double);									///<sets the range of the interaction potential
			void setNeighborStyle(int);											///<sets the neighbor style, or the way the cell neighbor list is built

			void clear();														///<Clears the contents of all the cells, their memory space, and the list of cells
			virtual void clearContents();										///<Clears the contents of all the cells
			void placeInCorrectCell(int, coord_t&);								///<If object has changed position, removes it from old cell, places it in new cell
			void insert(int, coord_t&);											///<Inserts an object into the cell list
			void remove(int);													///<Removes an object from the cell list
			
			bool good() const;													///<Indicates whether the cell list is valid
			bool fail() const;													///<Indicates whether the cell list is invalid
			
			bool end() const;
			void resetIterator();
			std::pair<int, int> nextPair();

			
			/**
			@brief Special function to get the neighbors of an object
			*/
			void getNeighborsOf(int, std::vector<int>&);

			/**
			@brief List of modes for setNeighborStyle
			
			For an MD-style cell list where interactions are pairwise, HALF is 
			typically used, since we only need half of the cell neighbors
			to compute the total force.  For an MC-style cell list, FULL is 
			typically used, since we often need to look at an individual 
			particle and compute the force over paritcles in all neighboring 
			cells.
			*/
			enum {FULL=0, HALF=1};

		protected:
			virtual void rebuild();												///<Builds the cell list	
			virtual void divideBoxIntoCells();									///<Divides the box into cells
			double _interaction_range;											///<The range of the interaction potential			
			double _min_cell_size;												///<The smallest possible cell dimension (usually _interaction_range).
			int _neighbor_style;												///<An integer indicating how to build cell neighbor list
			bool _force_rebuild;												///<If true, then rebuild the list no matter what next time Rebuild is called
			bool _is_good;														///<<Indicates whether the cell list is valid
			
			std::vector<Cell> _cell;											///<List of all cells in the cell list	
			std::map<int, int> _cell_containing;								///<The cell containing an object
			
			int _ncellyz;														///<number of cells in yz plane
			int _ncellx;														///<number of cells in x dim
			int _ncelly;														///<number of cells in y dim
			int _ncellz;														///<number of cells in z dim
			double _cell_size[3];												///<the size of a cell
			double _inv_cell_size[3];											///<inverse of the size of a cell 
			coord_t _box_min, _box_max;											///<The box boundary
			std::vector<bool> _periodic;										///<periodicity
		
		private:
			int _i, _j, _k, _n, _c;
			std::pair<int, int> _current_pair;
			bool _same_cell, _end;
			std::vector<int> _tag;
			std::map<int, int> _index_of_tag;									///<The index of a given tag in the _tag array
			
	};
}
	
#endif

