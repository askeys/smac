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
 *  StandardCell.cpp
 *  glotzilla
 *
 *  Created by Aaron Keys on 11/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "Cell.h"

namespace smac
{
	/**
	*/
	Cell::Cell()
	{
	}
	
	/**
	*/
	Cell::~Cell()
	{
	}

	/**
	clears the contents of the cell and the cell neighbors.
	*/
	void Cell::clear()
	{	
		_neighbor.clear();
		_tag.clear();
		_index_of.clear();		
	}

	/**
	clears the contents of the cells only (i.e., all tags contained in the cell).
	*/
	void Cell::clearContents()
	{		
		_tag.clear();
		_index_of.clear();		
	}		

#pragma mark GETS

	/**
	@return the number of tags in the cell
	*/
	int Cell::getNumberOfNeighbors() const
	{
		return _neighbor.size();
	}

	/**
	@return the number of tags in the cell
	*/
	int Cell::getNumberOfObjects()
	{
		return _tag.size();
	}

	/**
	*/
	std::vector<int> Cell::getNeighbors()
	{
		return _neighbor;
	}
	
	/**
	*/
	std::vector<int>  Cell::getContents()
	{
		return _tag;
	}

#pragma mark SETS	
#pragma mark OTHER
	
	/**
	@param tag is the tag to be inserted
	*/
	void Cell::insert(int tag)
	{				
		//save the index for fast removal
		_index_of[tag] = _tag.size();							
		_tag.push_back(tag);
	}

	/**
	@param tag is the tag to be removed
	*/
	void Cell::remove(int tag)
	{		
		int index_of_removed_tag = _index_of[tag];								//save the index of the tag to be removed
		int tmp_tag = _tag[_tag.size()-1];										//grab the last tag in the cell
		_tag[index_of_removed_tag] = tmp_tag;									//move that tag to the position of the removed tag
		_tag.pop_back();														//delete the last element, which is now a duplicate
		
		_index_of[tmp_tag] = index_of_removed_tag;								//save the index of the swapped tag
		_index_of.erase(tag);													//erase this key from the map		
	}	
		
#pragma mark NEIGHBOR OPERATIONS
	
	/**
	@param nbr is the index of a neighboring cell 
	*/
	void Cell::addNeighbor(int nbr)
	{
		_neighbor.push_back(nbr);
	}
}
