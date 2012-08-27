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
 *  MoveVelocityVerlet.h
 *  glotzilla
 *
 *  Created by Aaron Keys on 7/13/08.
 *  Copyright 2008 __MyCompanynatomame__. All rights reserved.
 *
 */

#include "cluster.h"
#include "space.h"

#include "CellList.h"

#include <iostream>
#include <cmath>
#include <map>
#include <vector>


#include <assert.h>
#include <cstdlib>

namespace smac
{
	/**
	\brief a rule to cluster points within a set range with per type cutuff
	\ingroup cluster
	*/
	bool clstrrulerangetype(int i, int j, arg_t ptr)
	{
		assert(ptr != 0x0 && ptr != NULL);
		clstrrulerangetype_info* args = (clstrrulerangetype_info*)ptr;

		if (args->box.periodic.size() == 0) {
			args->box.periodic.push_back(true);
			args->box.periodic.push_back(true);
			args->box.periodic.push_back(true);
		}
		if (args->box.period.size() == 0) {
			args->box.period.push_back(0.0);
			args->box.period.push_back(0.0);
			args->box.period.push_back(0.0);
		}
        if (args->rcutsq.size() == 0) {
            std::cerr << "***ERROR: clstrrulerangetype: per type cutoff";
            std::cerr << " (rcutsq) not defined." << std::endl;
            exit(1);
        }

		double rsq = distancesq(
			args->x[i], args->x[j], args->box.period, args->box.periodic);

		if (rsq < args->rcutsq[args->type[i]][args->type[j]]) {
			return true;
		}
		else {
			return false;
		}
	}


	/**
	\brief a rule to cluster points within a set range
	\ingroup cluster
	*/
	bool clstrrulerange(int i, int j, arg_t ptr)
	{
		assert(ptr != 0x0 && ptr != NULL);
		clstrrulerange_info* args = (clstrrulerange_info*)ptr;

		if (args->box.periodic.size() == 0) {
			args->box.periodic.push_back(true);
			args->box.periodic.push_back(true);
			args->box.periodic.push_back(true);
		}
		if (args->box.period.size() == 0) {
			args->box.period.push_back(0.0);
			args->box.period.push_back(0.0);
			args->box.period.push_back(0.0);
		}

		double rsq = distancesq(
			args->x[i], args->x[j], args->box.period, args->box.periodic);

		if (rsq < args->rcutsq) {
			return true;
		}
		else {
			return false;
		}
	}

	/**
	\brief a rule to cluster points within a shell
	\ingroup cluster
	*/
	bool clstrruleshell(int i, int j, arg_t ptr)
	{
		assert(ptr != 0x0 && ptr != NULL);
		clstrruleshell_info* args = (clstrruleshell_info*)ptr;

		if (args->box.periodic.size() == 0) {
			args->box.periodic.push_back(true);
			args->box.periodic.push_back(true);
			args->box.periodic.push_back(true);
		}
		if (args->box.period.size() == 0) {
			args->box.period.push_back(0.0);
			args->box.period.push_back(0.0);
			args->box.period.push_back(0.0);
		}

		double rsq = distancesq(
			args->x[i], args->x[j], args->box.period, args->box.periodic);

		if (rsq < args->rcutsqouter) {
			if (rsq >= args->rcutsqinner) {
				return true;
			}
			else {
				return false;
			}
		}
		else {
			return false;
		}
	}

	void add2cluster( int i, int j, clstrlist_t& cluster)
	{
		//if the cluster already exists, we can find it by looking it
		//up by the tag of its center atom
		if (cluster[i].size() == 0) {
			cluster[i].push_back(i);
		}
		if (cluster[j].size() == 0) {
			cluster[j].push_back(j);
		}
		cluster[i].push_back(j);
		cluster[j].push_back(i);
	}

	/**
	\brief Cluster points using a cell list for increased efficiency
	\ingroup cluster

	This algorithm loops over all pairs of points within a specified distance
	of one another and adds points to the cluster if they satisfy the clustering
	rule.  This increases computational efficientcy for short-range
	clustering rules.
	*/
	void clstrshortrange(coordlist_t& x, box_t& box, double range,
		clstrrule_t cluster_function, arg_t args, clstrlist_t& cluster)
	{
		int n = x.size();
		assert(n>0);

		cluster.clear();
		cluster.resize(n);

		if (box.boxlo.size() == 0) {
			box.boxlo.push_back(0.0);
			box.boxlo.push_back(0.0);
			box.boxlo.push_back(0.0);
		}
		if (box.boxhi.size() == 0) {
			box.boxhi.push_back(0.0);
			box.boxhi.push_back(0.0);
			box.boxhi.push_back(0.0);
		}

		CellList cell_list;
		cell_list.clear();
		cell_list.setBox(box.boxlo, box.boxhi);
		cell_list.setInteractionRange(range);
		cell_list.setNeighborStyle(CellList::HALF);

		if (cell_list.good()) {
			for(int i=0; i<n; i++) {
				cell_list.insert(i, x[i]);
			}
			cell_list.resetIterator();
			for (;;) {
				std::pair<int, int> p = cell_list.nextPair();
				if (cell_list.end()) {
					break;
				}
				if (cluster_function(p.first, p.second, args)) {
					add2cluster(p.first, p.second, cluster);
				}
			}
		}
		else {
			//std::cerr << "WARNING...\n";
			clstr(n, cluster_function, args, cluster);
		}
	}

	/**
	\brief a rule to cluster points within a set range
	\ingroup cluster

	This algorithm loops over all pairs of points adds points to the cluster if
	they satisfy the clustering rule
	*/
	void clstr(int n, clstrrule_t cluster_function, arg_t args, clstrlist_t& cluster)
	{
		assert(n > 0);
		cluster.clear();
		cluster.resize(n);

		for (int i=0; i<n; i++) {
			for (int j=i+1; j<n; j++) {
				if (cluster_function(i, j, args)) {
					add2cluster(i, j, cluster);
				}
			}
		}
	}

	void add2cluster_recursive(int i, int n, std::vector<int>& cluster,
		std::vector<bool>& counted,
		clstrrule_t cluster_function, arg_t args)
	{
		counted[i] = true;
		cluster.push_back(i);

		for (int j=0; j<n; j++) {
			if(!(counted[j])) {
				if (cluster_function(i, j, args)) {
					add2cluster_recursive(
						j, n, cluster, counted, cluster_function, args);
				}
			}
		}
	}

	void clstrhierarchical(int n, clstrrule_t cluster_function, arg_t args, clstrlist_t& cluster)
	{
		assert(n > 0);
		cluster.clear();

		std::vector<bool> counted(n, false);

		for(int i=0; i<n; i++) {
			if(!counted[i]) {
				std::vector<int> cluster_i;
				add2cluster_recursive(
					i, n, cluster_i, counted, cluster_function, args);
				cluster.push_back(cluster_i);
			}
		}
	}


	void add2cluster_recursive(int i, int n, std::vector<int>& cluster,
		std::vector<bool>& counted, CellList& cell_list,
		clstrrule_t cluster_function, arg_t args)
	{
		counted[i] = true;
		cluster.push_back(i);
		std::vector<int> nbr;
		cell_list.getNeighborsOf(i, nbr);

		for (unsigned int j=0; j<nbr.size(); j++) {
			if(!(counted[nbr[j]])) {
				if (cluster_function(i, nbr[j], args)) {
					add2cluster_recursive(
						nbr[j], n, cluster, counted, cell_list, cluster_function, args);
				}
			}
		}
	}


	void clstrshortrangehierarchical(coordlist_t& x,
		box_t& box,
		double range,
		clstrrule_t cluster_function, arg_t args,
		clstrlist_t& cluster)
	{
		int n = x.size();
		cluster.clear();
		cluster.resize(n);

		if (box.boxlo.size() == 0) {
			box.boxlo.push_back(0.0);
			box.boxlo.push_back(0.0);
			box.boxlo.push_back(0.0);
		}
		if (box.boxhi.size() == 0) {
			box.boxhi.push_back(0.0);
			box.boxhi.push_back(0.0);
			box.boxhi.push_back(0.0);
		}
            
		CellList cell_list;
		cell_list.clear();
		cell_list.setBox(box.boxlo, box.boxhi);
		cell_list.setInteractionRange(range);
		cell_list.setNeighborStyle(CellList::FULL);

		cluster.clear();
		std::vector<bool> counted(n, false);

		if (cell_list.good()) {
			for(int i=0; i<n; i++) {
				cell_list.insert(i, x[i]);
			}
			for(int i=0; i<n; i++) {
				if(!counted[i]) {
					std::vector<int> cluster_i;
					add2cluster_recursive(i, n, cluster_i,
						counted, cell_list, cluster_function, args);
					cluster.push_back(cluster_i);
				}
			}
		}
		else {
			//std::cerr << "WARNING...\n";
			clstrhierarchical(n, cluster_function, args, cluster);
		}
	}
}
