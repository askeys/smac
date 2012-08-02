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
 *  analysis.cpp
 *  smac
 *
 *  Created by askeys on 12/1/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "stdanalysis.h"
#include "cluster.h"
#include "preprocess.h"
#include "space.h"
#include <iostream>
#include "HistogramConst.hpp"
#include "io.h"
#include <stdexcept>


#include <assert.h>

namespace smac
{
	/**
	\brief computes a local descriptor for each particle
    \param x is a list of coordinates
    \param box defines the boundaries
    \param range defines the interaction range for what is "local"
    \param descfun is the shape descriptor to use
    \param descargs is a pointer to an argument struct for the descriptor
    \param cluster is a placeholder for the local clusters created for each
        particle
    \param descriptors is a placeholder for the list of output descriptors
    \param info is a struct containing additional optional information to pass
        to localorder

    This function creates a small local cluster about each point and computes
    a shape descriptor.  The cluster is translated
    */
	void localorder(coordlist_t& x, box_t& box, double range,
		shpdescfun_t descfun, arg_t descargs, clstrlist_t& cluster,
		shpdesclist_t& descriptors, localorder_info& info)
	{
		assert(x.size() > 0);
		int dim = x[0].size();
		descriptors.clear();

		clstrrulerange_info args;
		args.x = x;
		args.box = box;
		args.dim = dim;
		args.rcutsq = range*range;

		cluster.clear();
		clstrshortrange(x, box, range, clstrrulerange, &args, cluster);

		for (unsigned int i=0; i<cluster.size(); i++) {
			coordlist_t x_cluster(cluster[i].size());

			if (cluster[i].size() == 0) {
				cluster[i].push_back(i);
                shpdesc_t sd;
                sd.push_back(component_t(0.0));
                descriptors.push_back(sd);
                continue;
			}

			for (unsigned int j=0; j<cluster[i].size(); j++) {
				x_cluster[j] = x[cluster[i][j]];
			}
            if (info.origin == CENTROID) {
                unmapcentroid(x_cluster, box);
                transcentroid(x_cluster);
            }
            else if (info.origin == CENTER_PARTICLE) {
                coord_t x_c = x[cluster[i][0]];
                unmap(x_cluster, x_c, box);
            }
            else if (info.origin == SKIP_TRANSLATE) {
                //do nothing
            }
            else {
                std::cerr << "***ERROR: localorder: bad \"origin\" argument in";
                std::cerr << " localorder_info struct.\n";
                throw std::runtime_error("Error calling localorder");
            }
            //std::cerr << x_cluster << std::endl;
			shpdesc_t sd;
			weight_t f(x_cluster.size(), 1.0);
			shapedata_t shape_cluster(x_cluster, f);
			descfun(shape_cluster, sd, descargs);
			descriptors.push_back(sd);
		}
	}

	void localorder_central_atom(coordlist_t& x, box_t& box, double range,
		shpdescfun_t descfun, arg_t descargs, clstrlist_t& cluster,
		shpdesclist_t& descriptors)
	{
		assert(x.size() > 0);
		int dim = x[0].size();
		descriptors.clear();

		clstrrulerange_info args;
		args.x = x;
		args.box = box;
		args.dim = dim;
		args.rcutsq = range*range;

		cluster.clear();
		clstrshortrange(x, box, range, clstrrulerange, &args, cluster);

		for (unsigned int i=0; i<cluster.size(); i++) {
			coordlist_t x_cluster(cluster[i].size());
			if (cluster[i].size() == 0) {
				continue;
			}

			coord_t x_c;
			x_c.push_back(x[cluster[i][0]][0]);
			x_c.push_back(x[cluster[i][0]][1]);
			x_c.push_back(x[cluster[i][0]][2]);


			for (unsigned int j=0; j<cluster[i].size(); j++) {
				x_cluster[j] = x[cluster[i][j]];
			}

			unmap(x_cluster, x_c, box);
			//translate(x_cluster, x_c);
			shpdesc_t sd;
			weight_t f(x_cluster.size(), 1.0);
			shapedata_t shape_cluster(x_cluster, f);

			descfun(shape_cluster, sd, descargs);
			descriptors.push_back(sd);
		}
	}

	void bod_central_atom(coordlist_t& x,  box_t& box, double range,
		shpdescfun_t descfun, arg_t descargs, clstrlist_t& cluster,
		shpdesc_t& descriptor, coordlist_t& x_combined)
	{
		assert(x.size() > 0);
		int dim = x[0].size();
		descriptor.clear();

		//coordlist_t x_combined;

		clstrrulerange_info args;
		args.x = x;
		args.box = box;
		args.dim = dim;
		args.rcutsq = range*range;

		cluster.clear();
		clstrshortrange(x, box, range, clstrrulerange, &args, cluster);


		for (unsigned int i=0; i<cluster.size(); i++) {
			coordlist_t x_cluster(cluster[i].size());
			if (cluster[i].size() == 0) {
				continue;
			}

			coord_t x_c;
			x_c.push_back(x[cluster[i][0]][0]);
			x_c.push_back(x[cluster[i][0]][1]);
			x_c.push_back(x[cluster[i][0]][2]);


			for (unsigned int j=0; j<cluster[i].size(); j++) {
				x_cluster[j] = x[cluster[i][j]];
			}

			unmap(x_cluster, x_c, box);
			for (unsigned int j=1; j<cluster[i].size(); j++) {
				x_combined.push_back(x_cluster[j]);
			}


		}
		weight_t f(x_combined.size(), 1.0);
		shapedata_t shape_cluster(x_combined, f);
		descfun(shape_cluster, descriptor, descargs);
	}


	/**
	\brief determines if each particle is in a liquid or crystal configuration
	\ingroup xstal
	\ingroup analysis
	\param x is the coordinates
	\param box is the box
	\param rangeop is the order parameter range for a coordination shell
	\param rangematch is the correlation range for comparing coordination shells
	\param bondcut is the threshold for what is considered a "solid bond"
	\param fbondcut is the minimum number of solid bonds for a solid particle
	\note rangeop is typically equal to rangematch

	The main idea here is that particles in a crystal have coodination shells
	that are correlated with those of their neighbors in terms of shape and
	orientation.  If the correlation between the order parameters for a particle
	and its neighbor exceed a certain value, we consider this a "solid bond."
	Since particles might randomly have solid-like bonds with their neighbors
	in the liquid, we only consider particles that exceed a certain fraction of
	solid bonds as being solid-like.
	*/
	void crystalorliquid(coordlist_t& x, box_t& box,
		double rangeop, double rangematch,
		shpdescfun_t descfun, arg_t descargs, matchfun_t matchfun,
		double bondcut, double fbondcut, std::vector<bool>& solid)
	{
		assert(x.size() > 0);

		shpdesclist_t sd;
		clstrlist_t cluster;
		solid.resize(x.size());
        localorder_info loinfo;
		localorder(x, box, rangeop, descfun, descargs, cluster, sd, loinfo);
		double rangematchsq = rangematch*rangematch;

		if (rangematch > rangeop) {
			clstrrulerange_info args;
			args.x = x;
			args.box = box;
			args.dim =  x[0].size();
			args.rcutsq = rangematch*rangematch;

			cluster.clear();
			clstrshortrange(x, box, rangematch, clstrrulerange, &args, cluster);
		}

		for (unsigned int i=0; i<cluster.size(); i++) {
			int nbond = 0;
			for (unsigned int j=1; j<cluster[i].size(); j++) {
				int neighbor = cluster[i][j];
				if (distancesq(x[i], x[neighbor], box.period, box.periodic)
					< rangematchsq) {
					if (sd[i].size() != sd[neighbor].size()) {
						continue;
					}
					double m = matchfun(sd[i], sd[neighbor], 0.0, 1.0);
					if (m > bondcut) {
						nbond++;
					}
				}
			}
			solid[i] = false;
			if (cluster[i].size() > 1) {
				if (nbond/(double)(cluster[i].size()-1) >= fbondcut) {
					solid[i] = true;
				}
			}
		}
	}

	void crystalgrains(coordlist_t& x, box_t& box,
		crystalanalysis_info& ca, clstrlist_t& grains)
	{
		assert(x.size() > 0);

		double rangeop = ca.rangedescriptor;
		double rangematch = ca.rangecompare;
		double rangeclstr = ca.rangecluster;
		shpdescfun_t descfun = ca.descriptor;
		arg_t descargs = ca.descriptorarg;
		matchfun_t matchfun = ca.matchfun;
		double bondcut = ca.crystalbondcut;
		double fbondcut = ca.crystalparticlecut;		
        
		//save the descriptors if necessary, otherwise don't
		shpdesclist_t sd1, *sdp;
		if (ca.savedescriptors) {
			sdp = &(ca.saveddescriptors);
		}
		else {
			sdp = &sd1;
		}
		shpdesclist_t& sd = *sdp;

		//save the clusters if necessary, otherwise don't
		clstrlist_t cluster1, *clusterp;
		if (ca.saveclusters) {
			clusterp = &(ca.savedclusters);
		}
		else {
			clusterp = &cluster1;
		}
		clstrlist_t& cluster = *clusterp;

		//done initializing stuff
		coordlist_t x_crystal;
		std::vector<int> index_crystal;
        localorder_info loinfo;
        loinfo.origin = CENTER_PARTICLE;
		localorder(x, box, rangeop, descfun, descargs, cluster, sd, loinfo);					
		double rangematchsq = rangematch*rangematch;
        
		if (rangematch > rangeop) {
			clstrrulerange_info args;
			args.x = x;
			args.box = box;
			args.dim =  x[0].size();
			args.rcutsq = rangematch*rangematch;

			cluster.clear();
			clstrshortrange(x, box, rangematch, clstrrulerange, &args, cluster);
		}

		for (unsigned int i=0; i<cluster.size(); i++) {
			int nbond = 0;

			for (unsigned int j=1; j<cluster[i].size(); j++) {
				int neighbor = cluster[i][j];
				if (distancesq(x[i], x[neighbor], box.period, box.periodic)
					< rangematchsq) {
					if (sd[i].size() != sd[neighbor].size()) {
						continue;
					}
                    double m = matchfun(sd[i], sd[neighbor], -1.0, 1.0);
					if (ca.matchmode == MATCHNEGATE) {
                        m *= -1.0;
                    }
					if (ca.matchmode == MATCHABS) {
                        m = fabs(m);
                    }
                    if (m > bondcut) {
						nbond++;
					}
				}
			}
			if (cluster[i].size() > 1) {
                if (fbondcut <= 1.0) {
                    if (nbond/(double)(cluster[i].size()-1) >= fbondcut) {
                        x_crystal.push_back(x[i]);
                        index_crystal.push_back(i);
                    }
                }
                else {
                    if (nbond >= fbondcut) {
                        x_crystal.push_back(x[i]);
                        index_crystal.push_back(i);                        
                    }
                }
			}
		}

		clstrlist_t temp_grains;
		clstrrulerange_info cargs;
		cargs.box = box;
		cargs.rcutsq = rangeclstr*rangeclstr;
		cargs.x = x_crystal;
		clstrshortrangehierarchical(
			x_crystal, box, rangeclstr, clstrrulerange, &cargs, temp_grains);

		grains.resize(temp_grains.size());
		for (unsigned int i=0; i<temp_grains.size(); i++) {
			clstr_t grain(temp_grains[i].size());
			for (unsigned int j=0; j<temp_grains[i].size(); j++) {
				grain[j] = index_crystal[temp_grains[i][j]];
			}
			grains[i] = grain;
		}
	}

	void testcrystalcriterion(coordlist_t& x, box_t& box,
		double rangeop, double rangematch, shpdescfun_t descfun,
		arg_t descargs, matchfun_t matchfun,
		double bondcut, function1d_t& dothistf, function1d_t& bondhistf,
        double histstep)
	{
		shpdesclist_t sd;
		clstrlist_t cluster;

        localorder_info loinfo;
		localorder(x, box, rangeop, descfun, descargs, cluster, sd, loinfo);
		HistogramConst<double> bondhist(0, histstep, 1.0);
		HistogramConst<double> dothist(0, histstep, 1.0);
		double rangematchsq = rangematch*rangematch;

		for (unsigned int i=0; i<cluster.size(); i++) {
			int nbond = 0;
			for (unsigned int j=1; j<cluster[i].size(); j++) {
				int neighbor = cluster[i][j];
				if (distancesq(x[i], x[neighbor], box.period, box.periodic)
					< rangematchsq) {
					if (sd[i].size() != sd[neighbor].size()) {
						continue;
					}
					double m = matchfun(sd[i], sd[neighbor], 0.0, 1.0);
					dothist.insert(m);
					if (m > bondcut) {
						nbond++;
					}
				}
			}
			if (cluster[i].size() > 1) {
				bondhist.insert(nbond/(double)(cluster[i].size()-1));
			}
		}

		int nbin = bondhist.getNumberOfBins();
		dothistf.clear();
		bondhistf.clear();
		std::pair<double, double> bin;

		for (int i=0; i<nbin; i++) {
			bin = dothist.getBin(i);
			dothistf[bin.first] = bin.second;
			bin = bondhist.getBin(i);
			bondhistf[bin.first] = bin.second;
		}
	}
}
