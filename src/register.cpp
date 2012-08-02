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
 *  register.cpp
 *  smac
 *
 *  Created by askeys on 12/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "register.h"
#include "preprocess.h"
#include "space.h"
#include "rand.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#ifdef _USE_ICP
#include "../3rdparty/icp/icpCpp.h"
#endif
#ifdef _USE_PCA
#include "../3rdparty/pca/pca.h"
#endif

#include "shpdescrms.h"
#include "matchfun.h"

namespace smac
{	

    void icptrialuniform(coordlist_t& x, int iter, int maxiter)
    {
        int rpd = pow(maxiter, 1.0/3.0);    //rotations per dimension
        
        int iterz = iter%rpd;
        double rotz = iterz/(double)rpd*2.0*M_PI;
        int itery = iter/rpd;
        double roty = itery/(double)rpd*2.0*M_PI;        
        int iterx = iter/(rpd*rpd);
        double rotx = iterx/(double)rpd*2.0*M_PI;
        rotate(x, rotx, roty, rotz);
    }

    /**
    
    */
    void icptrialrandom(coordlist_t& x)
    {
        std::vector<double> q(4);
        randquaternion(q);
        rotateq(x, q);
    }

#ifdef _USE_ICP
    void callicp(coordlist_t& x_query, coordlist_t& x_ref, 
        Tree* tree, std::vector<unsigned int> randvec, std::vector<double> weights,
        int niter, double tol)
    {    
		int n = x_query.size(); 
		int n_ref = x_ref.size();

        //1D arrays for ICP
		std::vector<double> mem_x(n*3);
		for (int i=0; i<n; i++) {
			mem_x[i*3+0] = x_query[i][0];
			mem_x[i*3+1] = x_query[i][1];
			mem_x[i*3+2] = x_query[i][2];
		}
		std::vector<double> mem_ref(n_ref*3);
		for (int i=0; i<n_ref; i++) {
			mem_ref[i*3+0] = x_ref[i][0];
			mem_ref[i*3+1] = x_ref[i][1];
			mem_ref[i*3+2] = x_ref[i][2];
		}
		
		int sizernd = ceil(1.45*n);
        //possibly change code to accept weights			
        double tr[9] = {0.0}, tt[3] = {0.0};

        icp(tr, tt, 
            &mem_ref[0], n_ref, 
            &mem_x[0], &weights[0], n, 
            &randvec[0], randvec.size(), sizernd,
            niter, tree);
                
        for (int i=0; i<n; i++) {
                
            double temp[3];
            temp[0] = tr[0*3+0] * x_query[i][0] + 
                tr[1*3+0] * x_query[i][1] + tr[2*3+0] * x_query[i][2];
            temp[1] = tr[0*3+1] * x_query[i][0] + 
                tr[1*3+1] * x_query[i][1] + tr[2*3+1] * x_query[i][2];
            temp[2] = tr[0*3+2] * x_query[i][0] + 
                tr[1*3+2] * x_query[i][1] + tr[2*3+2] * x_query[i][2];
        
            x_query[i][0] = temp[0];
            x_query[i][1] = temp[1];
            x_query[i][2] = temp[2];
        }
        
        for (int i=0; i<n; i++) {
            x_query[i][0] += tt[0]; 
            x_query[i][1] += tt[1]; 
            x_query[i][2] += tt[2];
        }
	}
#endif
    
    /**
    \brief register objects using iterative-closest-point (ICP)
    \ingroup reg
    \param x_query is the coordinates of the object to be registered (will be 
        modified)
    \param x_ref is the coordinates of the object to register with (will not be 
        modified)
    \param arg is a pointer to a registericp_info struct containing additional 
        optional arguments
    */
	void registericp(coordlist_t& x_query, coordlist_t& x_ref, arg_t arg)
	{
#ifdef _USE_ICP
        registericp_info* info = static_cast<registericp_info*> (arg);
        registericp_info info_tmp;
        if (info == 0x0) {
            info = &info_tmp;
        }
        
		int n = x_query.size(); 
		int n_ref = x_ref.size();

        //some variables to only initialize once:
        bool icp_initialized = false;
        Tree* tree = 0x0; 
		std::vector<double> weights(n, 1.0);
		std::vector<unsigned int> randvec(n);
        rmsassign_info assignargs;
        rmsdesc_info rmsargs;
        
        shpdesc_t descriptor_ref, descriptor_query;
        shapedata_t shape_query(x_query);
        shapedata_t shape_ref(x_ref);
        rmsdesc(shape_ref, descriptor_ref, &rmsargs);

        double max_match = -1e20;
        coordlist_t x_best = x_query;
        coordlist_t x_original = x_query;

        for (int iter = 0; iter<=info->maxrotations; iter++) {
            //the first time, leave it as-is, afterward try rotating the object
            if (iter > 0) {
                x_query = x_original;
                if (info->trialtype == RANDOM) {
                    icptrialrandom(x_query);
                }
                else if (info->trialtype == LINEAR) {
                    icptrialrandom(x_query);
                }
                else {
                    std::cerr << "***ERROR: icp rotation trialtype undefined\n";
                    throw std::runtime_error("ERROR using ICP");
                }
            }
        
            shape_query.x = x_query;
            rmsdesc(shape_query, descriptor_query, &rmsargs);
            rmsfastassign(descriptor_query, descriptor_ref, &assignargs);

            double match = info->matchfun(descriptor_query, descriptor_ref, 0.0, 1.0);
            if (match > 1.0-info->tol) {
                return;
            }
            if (match > max_match) {
                x_best = x_query;
                max_match = match;                        
			}
            
            //dont even bother trying icp if we can't even get somewhat close
            if (match > info->icptol) {
                if (!icp_initialized) {
                    std::vector<double> x_ref_mem(n_ref*3);
                    //sideways arrays for kdtree (annoying)
                    int xoffset=0, yoffset=n_ref, zoffset=n_ref*2;
                    for (int i=0; i<n_ref; i++) {
                        x_ref_mem[i+xoffset] = x_ref[i][0];
                        x_ref_mem[i+yoffset] = x_ref[i][1];
                        x_ref_mem[i+zoffset] = x_ref[i][2];
                    }

                    //index array
                    std::vector<int> index(n_ref);
                    for (unsigned int i=0; i<index.size(); i++) {
                        index[i] = i;
                    }
                            
                    tree = build_kdtree(&x_ref_mem[0], 
                        n_ref, 3, &index[0], n_ref, 0);
                    tree->dims = 3;
                    //display_tree(tree->rootptr, 3);
                    
                    std::vector<bool> counted(n, false);
                    int count = 0;
                    while (count < n) {
                        int t = rand() % n;
                        if (!counted[t]) {
                            randvec[count] = t;
                            counted[t] = true;
                            count++;
                        }
                    }
                    icp_initialized = true;
                }                
                callicp(x_query, x_ref, tree, randvec, weights, info->maxiter, info->tol);
            }
            else {
                continue;
            }
            //sometimes we get NAN
            if (x_query[0][0] != x_query[0][0]) {
                continue;
            }
            shape_query.x = x_query;    
            rmsdesc(shape_query, descriptor_query, &rmsargs);
            rmsfastassign(descriptor_query, descriptor_ref, &assignargs);
            match = info->matchfun(descriptor_query, descriptor_ref, 0.0, 1.0);
            if (match > 1.0-info->tol) {
                return;
            }
            if (match > max_match) {
                x_best = x_query;
                max_match = match;                        
			}
		}
        if (icp_initialized) {
            free_tree(tree->rootptr);
        }
        x_query = x_best;
#else
        std::cerr << "***ERROR: attempt to use ICP registration without ";
        std::cerr << "enabling the ICP library." ;
        throw std::runtime_error("Error using ICP");
#endif
    }






    /**
    \brief register objects using principle components analysis (PCA)
    \ingroup reg
    \param x is the coordinates of the object
    \param arg is a pointer to a registerpca_info struct containing additional 
        optional arguments
    */
	void registerpca(coordlist_t& x, arg_t arg)
	{
#ifdef _USE_PCA
        registerpca_info* info = static_cast<registerpca_info*> (arg);
        registerpca_info info_tmp;
        if (info == 0x0) {
            info = &info_tmp;
        }
        char mode = info->mode;
        
		int n = x.size();
		if (n == 0) {
			return;
		}
		int m = x[0].size();
		
		std::vector<float> x1d(n*m);
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				x1d[i*m+j] = x[i][j];
			}
		}
		
		pca(&x1d[0], n, m, mode);
		
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				x[i][j] = x1d[i*m+j];
			}
		}		
#else
        std::cerr << "***ERROR: attempt to use PCA registration without ";
        std::cerr << "enabling the PCA library." ;
        throw std::runtime_error("Error using PCA");
    
#endif
	}
    
    /**
    \brief Aligns the object with the real component of the axis of max symmetry
    \ingroup reg
    \param x_orig is the input cluster to be registered
    \param info is a fourier_info struct containing parameters for the fourier
        descriptor
    */
    double registerzerocomplexcomonent(coordlist_t& x_orig, fourier_info& info)
    {
        coordlist_t x = x_orig;
        for (unsigned int i=0; i<x.size(); i++) {
            x[i][2] = 0.0;
        }
        clstrlist_t clusters;
        for (unsigned int j=0; j<x.size(); j++) {
            for (int k=0; k<3; k++) {
                x[j][k] -= x[0][k];
            }
        }
        shpdesc_t sd;
        shapedata_t shape(x);
        fourierdesc2d(shape, sd, &info);
        double xt = sd[0].real();
        double yt = sd[0].imag();
        double ltheta = atan2(yt, xt);
        return (-ltheta / (double)info.frequency[0]);
    }
    
    /**
    \brief registers an object with the axis of max symmetry
    \ingroup reg
    \param x_orig is the input cluster to be registered
    \param q_best is a placeholder for the quaternion that maximizes symmetry
    \param info is a registersymmetryaxis_info struct containing parameters for
        the function.
    This routine is used for registering 3D crystals in the xy plane.  This is 
    useful for calculating the diffraction image or some other plane-dependent
    descriptor.  The current routine uses simulated annealing to converge 
    to a local minimum; thus the result may not be the true axis of max 
    symmetry.  However, this is not necessarily problematic, since the true
    axis is often a trivial (often 90 degree) rotation of the other axis.  The
    routine does not actually register the coordinates; rather it returns the 
    quaternion that maximizes the symmetry.  The structure can be rotated using
    the rotateq function.  The function is set up this way so that the user can
    use a small cluster from the system to find the best rotation (usually only
    ~100 particles are needed, since all the particles in a cystal are aligned) 
    and then use the quaternion to orient the whole crystal (which may consist 
    of up to millions of particles).
    */
    void registersymmetryaxis(
        coordlist_t& x_orig, 
        std::vector<double>& q_best, 
        registersymmetryaxis_info& info)
    {        
        if (info.fourierargs == 0x0) {
            std::cerr << "***ERROR: registersymmetryaxis: ";
            std::cerr << "registersymmetryaxis_info member \"fourierargs\"";
            std::cerr << "must be defined.";
            throw std::runtime_error("Error using registersymmetryaxis");
        }
        
        q_best = std::vector<double>(4, 0.5);
        double e_best = 1e100;
        
        for (int s=0; s<info.nstartpoints; s++) {
            std::vector<double> q(4, 0.0);
            randquaternion(q);
            double e_best_i = 0.0;
            std::vector<double> q_best_i = q;
            double dqmax = info.dqmax;
            for (int i=0; i<info.niter; i++) {
                std::vector<double> qr(4, 0.0);
                q = q_best_i;
                randquaternion(qr);
                q[0] += dqmax*qr[0];
                q[1] += dqmax*qr[1];
                q[2] += dqmax*qr[2];
                q[3] += dqmax*qr[3];
                double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
                q[0] /= norm;
                q[1] /= norm;
                q[2] /= norm;
                q[3] /= norm;
                
                coordlist_t x = x_orig;
                rotateq(x, q);
                                
                for (unsigned int j=0; j<x.size(); j++) {
                    x[j][2] = 0.0;
                }
                for (unsigned int j=0; j<x.size(); j++) {
                    for (int k=0; k<3; k++) {
                        x[j][k] -= x[0][k];
                    }
                }
                shpdesc_t sd;
                shapedata_t shape(x);
                fourierdesc2d(shape, sd, info.fourierargs);
                                                
                double energy = 0.0;
                for (unsigned int k=0; k<sd.size(); k++) {
                    energy -= abs(sd[k]);
                }
                if (exp((info.beta+i)*(e_best_i-energy)) >= drand48()) {
                    e_best_i = energy;
                    q_best_i = q;
                }
                else {
                    dqmax *= 0.95;
                }
                if (e_best_i > 0.5*e_best) {
                    break;
                }
            }
            if (e_best_i < e_best) {
                e_best = e_best_i;
                q_best = q_best_i;
            }
            if (e_best_i < info.eterminate) {
                break;
            }
        }    
    }
}