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
 *  io.cpp
 *  libsmac
 *
 *  Created by askeys on 12/4/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "io.h"
#include "space.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <map>
#include <assert.h>

namespace smac
{	

	/**
	\brief Load raw coordinates
	\param filename is the filename
	\param x is a list to hold the coordinates
	*/
	void load(const char* filename, coordlist_t& x, arg_t arg)
	{
		std::ifstream file(filename);
		if (file.fail()) {
			std::cerr << "Error: load: can't open file " 
				<< filename << ".\n";
			exit(1);
		}
		
		std::string str;
		while (std::getline(file, str)) {
			std::istringstream iss(str.c_str());
			double temp;
			coord_t xi;
			while (iss >> temp) {
				xi.push_back(temp);
			}
			if (xi.size() != 0) {
				x.push_back(xi);
			}
		}
	}
			
	void deltype(coordlist_t& x, std::vector<std::string>& types, std::string t)
	{
		assert(x.size() == types.size());

		coordlist_t xtmp = x;
		std::vector<std::string> typetmp = types;
				
		types.clear();
		x.clear();
		
		for (unsigned int i=0; i<xtmp.size(); i++) {
			if (typetmp[i] != t) {
				x.push_back(xtmp[i]);
				types.push_back(typetmp[i]);
			}
		}
	}

    /**
    \brief load data from a rasmol xyz file
    \ingroup io
    \param filename is the name of the xyz file
    \param xyz is a xyzfile struct to hold the data from the xyz file 
    */
	void loadxyz(const char* filename, xyzfile& xyz)
	{
		std::ifstream file(filename);
		assert(file.good());				
		loadxyz(file, xyz);
        file.close();
	}

    /**
    \brief load data from a rasmol xyz file
    \ingroup io
    \param is is a stream pointing to an xyz file
    \param xyz is a xyzfile struct to hold the data from the xyz file 
    */
	void loadxyz(std::istream& is, xyzfile& xyz)
    {
        int N=0;
        is >> N;
        xyz.n = N;
        std::getline(is, xyz.commentstr);
        std::getline(is, xyz.commentstr);
        xyz.x.resize(N, std::vector<double>(3,0.0));
        xyz.type.resize(N);
        for (int i=0; i<N; i++) {
            double x, y, z;
            std::string typestr;
            is >> typestr >> x >> y >> z;
            xyz.type[i] = typestr;
            xyz.x[i][0] = x;
            xyz.x[i][1] = y;
            xyz.x[i][2] = z;
        }
    }
	
	void loadxyz(const char* filename, coordlist_t& x, arg_t arg)
	{
		xyz_info* info = (xyz_info*)arg;
		xyz_info temp_info;
		if (arg == 0x0) {
			info = &temp_info;
		}
		
		std::istream* is;
		std::ifstream ifs;  
		if (info->instream != 0x0) {
			is = info->instream;
		}
		else {
			ifs.open(filename);
			is = &ifs;
		}
		
		assert(is->good());
		
		int n = 0;
		std::string str;
		std::getline(*is, str);
		std::istringstream iss(str.c_str());
		iss >> n;
		std::getline(*is, str);
		x.clear();
		info->type.clear();
		
		for (int i=0; i<n; i++) {
			std::string typei;
			double xi, yi, zi;
			*is >> typei >> xi >> yi >> zi;
			info->type.push_back(typei);
			coord_t xyz(3);
			xyz[0] = xi; xyz[1] = yi; xyz[2] = zi;
			x.push_back(xyz);
		}
		std::getline(*is, str);
		
		if (info->instream == 0x0) {
			ifs.close();
		}
	}

	/**
	\brief save an xyz file with variate types
	\ingroup io
	*/
	void savevarxyz(const char* filename, coordlist_t&x, xyz_info& info)
	{
		assert(info.type.size() > 0);
		assert(info.resevoir.size() > 0);
		std::ofstream xyzfile(filename);

		std::ostream *os;
		if (info.outstream == NULL) {
			assert(xyzfile.good());
			os = &xyzfile;
		}	
		else {
			assert(info.outstream->good());
			os = info.outstream;
		}

		int n=0;
		typedef std::map< std::string, int > resevoir_t;
		resevoir_t::iterator r;

		for (r = info.resevoir.begin(); r!=info.resevoir.end(); ++r) {
			n += r->second;
		}

		*os << n <<"\n\n";
				
		for (r = info.resevoir.begin(); r!=info.resevoir.end(); ++r) {
			std::string type = r->first;
			int npad = r->second;
			for (unsigned int i=0; i<x.size(); i++) {
				if (info.type[i] == type) {
					*os << type <<"\t" << x[i][0] <<"\t"<< x[i][1] <<"\t"<< x[i][2] <<"\n";
					npad--;
				}
			}
			for (int i=0; i<npad; i++) {
				*os << type <<"\t" << info.xres[0] <<"\t"<< info.xres[1] <<"\t"<< info.xres[2] <<"\n";			
			}
		}
		xyzfile.close();		
	}
		
	/**
	\brief save an xyz file
	\ingroup io
	*/
	void savexyz(const char* filename, coordlist_t&x, xyz_info& info)
	{	
		int n = x.size();

		std::ofstream xyzfile(filename);

		std::ostream *os;
		if (info.outstream == NULL) {
			assert(xyzfile.good());
			os = &xyzfile;
		}	
		else {
			assert(info.outstream->good());
			os = info.outstream;
		}

		*os << n <<"\n\n";
		if (info.type.size() > 0) {
			for (int i=0; i<n; i++) {
				*os << info.type[i] <<"\t" << x[i][0] <<"\t"<< x[i][1] <<"\t"<< x[i][2] <<"\n"; 
			}
		}
		else {
			for (int i=0; i<n; i++) {
				*os << "H\t" << x[i][0] <<"\t"<< x[i][1] <<"\t"<< x[i][2] <<"\n"; 
			}
		}
		xyzfile.close();			
	}



    /**
    \brief load data from a lammpstrj file
    \ingroup io
    \param filename is the name of the lammpstrj file
    \param lammpstrj is a lammpstrjfile struct to hold the data from the 
        lammpstrj file 
    */
	void loadlammpstrj(const char* filename, lammpstrjfile& lammpstrj)
	{
		std::ifstream file(filename);
		assert(file.good());				
		loadlammpstrj(file, lammpstrj);
        file.close();
	}

    /**
    \brief load data from a LAMMPS dump file
    \ingroup io
    \param is is a stream pointing to an LAMMPS dump file
    \param lammpstrj is a lammpstrjfile struct to hold the data from the 
        lammpstrj file
    \bug if velocities, bonds, etc. are in dump file, reading a long
        concatenated lammpstrj will fail
    */
	void loadlammpstrj(std::istream& is, lammpstrjfile& lammpstrj)
    {
		assert(is.good());
		
		int tempi, n;
		double t1, t2, t3, t4, t5, t6;
		char spc;
		std::string temp;
		//Timestep
		std::getline(is, temp);
		is >> tempi;
        lammpstrj.timestep = tempi;
		is >> spc;
		//Number of Atoms
		std::getline(is, temp);
		is >> n;
        lammpstrj.n = n;
		is >> spc;
		//Box
		std::getline(is, temp);
		is >> t1 >> t2 >> t3 >> t4 >> t5 >> t6;
		is >> spc;
		lammpstrj.box.boxlo.resize(3);
		lammpstrj.box.boxhi.resize(3);
		lammpstrj.box.period.resize(3);		
		lammpstrj.box.period[0] = t2-t1;
		lammpstrj.box.period[1] = t4-t3;
		lammpstrj.box.period[2] = t6-t5;
		lammpstrj.box.boxlo[0] = -0.5*(t2-t1);
		lammpstrj.box.boxlo[1] = -0.5*(t4-t3);
		lammpstrj.box.boxlo[2] = -0.5*(t6-t5);
		lammpstrj.box.boxhi[0] = 0.5*(t2-t1);
		lammpstrj.box.boxhi[1] = 0.5*(t4-t3);
		lammpstrj.box.boxhi[2] = 0.5*(t6-t5);

        if (is.fail()) {
            std::cerr << "***ERROR: bad file header for lammpstrj\n";
            throw std::runtime_error("Error reading lammpstrj\n");
        }                
        //Atoms		
		lammpstrj.type.resize(n);
		lammpstrj.x.resize(n);
        
		is >> spc;
        std::getline(is, temp);
		int particle_id;
				
        coord_t xi(3);
		
        for (int i=0; i<n; i++) {
			is >> particle_id >> tempi >> t1 >> t2 >> t3;
            lammpstrj.type[particle_id-1] = tempi;
			xi[0] = t1*lammpstrj.box.period[0] + lammpstrj.box.boxlo[0];
			pbc(xi[0], lammpstrj.box.period[0], lammpstrj.box.periodic[0]);
            xi[1] = t2*lammpstrj.box.period[1] + lammpstrj.box.boxlo[1];
			pbc(xi[1], lammpstrj.box.period[1], lammpstrj.box.periodic[1]);
			xi[2] = t3*lammpstrj.box.period[2] + lammpstrj.box.boxlo[2];
			pbc(xi[2], lammpstrj.box.period[2], lammpstrj.box.periodic[2]);
            lammpstrj.x[particle_id-1] = xi;
		}		
		is >> spc;
    }
    
	/**
	\brief Load a lammps dump file from a stream (typically a pre-opened file)
	\ingroup io
	\param x is the filename
	\param x is a list to hold the coordinates
	\param arg is a pointer to additional arguments (optional)
	\param box is a struct to hold the box information
	
	This function is typically used when we have a lammps trajectory file
	*/
	void loadlmptraj(const char* filename, coordlist_t& x, arg_t arg)
	{
		lmptraj_info* info = (lmptraj_info*)arg;
		
		std::istream *is;
		std::ifstream ifs;
		if (info->instream != 0x0) {
			is = info->instream;
		}
		else {
			ifs.open(filename);
			is = &ifs;
		}
		
		assert(is->good());
		
		int tempi, n;
		double t1, t2, t3, t4, t5, t6;
		char spc;
		std::string temp;
		//Timestep
		std::getline(*is, temp);
		*is >> tempi;
		*is >> spc;
		//Number of Atoms
		std::getline(*is, temp);
		*is >> n;
		*is >> spc;
		//Box
		std::getline(*is, temp);
		*is >> t1 >> t2 >> t3 >> t4 >> t5 >> t6;
		*is >> spc;
		info->box.boxlo.resize(3);
		info->box.boxhi.resize(3);
		info->box.period.resize(3);
		info->box.boxlo[0] = t1;
		info->box.boxlo[1] = t3;
		info->box.boxlo[2] = t5;
		info->box.boxhi[0] = t2;
		info->box.boxhi[1] = t4;
		info->box.boxhi[2] = t6;
		
		info->box.period[0] = t2-t1;
		info->box.period[1] = t4-t3;
		info->box.period[2] = t6-t5;
		//Atoms		
		info->type.clear();
		std::getline(*is, temp);
		int particle_id;
		
		coordlist_t x_temp(n);
		std::vector<std::string> type_temp(n);
		
		for (int i=0; i<n; i++) {
			*is >> particle_id >> tempi >> t1 >> t2 >> t3;
			std::ostringstream tmp;
			tmp << tempi;
			type_temp[particle_id-1]  = tmp.str();

			//info->type.push_back(tmp.str());
			coord_t xi;
			xi.push_back(t1*info->box.period[0] + info->box.boxlo[0]);
			xi.push_back(t2*info->box.period[1] + info->box.boxlo[1]);
			xi.push_back(t3*info->box.period[2] + info->box.boxlo[2]);
			//x.push_back(xi);
			for(int k=0; k<3; k++)
			{
				if(xi[k] > info->box.boxhi[k])
					xi[k] -= info->box.period[k];
				if(xi[k] < info->box.boxlo[k])
					xi[k] += info->box.period[k];
				
			}
			x_temp[particle_id-1] = xi;

		}
		
		for (int i=0; i<n; i++)
		{
			x.push_back(x_temp[i]);
			info->type.push_back(type_temp[i]);
		}
		
		*is >> spc;
		
		if (info->instream == 0x0) {
			ifs.close();
		}
		
	}
	
    /**
	\brief Load a lammps dump file from a stream (typically a pre-opened file) where we properly order the particles by id.
	\ingroup io
	\param is is the input stream to read
	\param x is a list to hold the coordinates
	\param type is a list to hold the atom types
	\param box is a struct to hold the box information
	
	This function is typically used when we have a lammps trajectory file. 
	Since order isn't guaranteed in the dump file, this function reorders by particle id.
	*/
	void loadlmptraj(
		std::istream& is, coordlist_t& x, arg_t arg)
	{
		lmptraj_info* info = (lmptraj_info*)arg;
	
		std::istream& lammps_data_file = is;
		assert(lammps_data_file.good());
		int tempi, n, particle_id;
		double t1, t2, t3, t4, t5, t6;
		char spc;
		std::string temp;
		//Timestep
		std::getline(is, temp);
		is >> tempi;
		is >> spc;
		//Number of Atoms
		std::getline(is, temp);
		is >> n;
		is >> spc;
		//Box
		std::getline(is, temp);
		is >> t1 >> t2 >> t3 >> t4 >> t5 >> t6;
		is >> spc;
		info->box.boxlo.resize(3);
		info->box.boxhi.resize(3);
		info->box.period.resize(3);
		info->box.boxlo[0] = t1;
		info->box.boxlo[1] = t3;
		info->box.boxlo[2] = t5;
		info->box.boxhi[0] = t2;
		info->box.boxhi[1] = t4;
		info->box.boxhi[2] = t6;
		
		info->box.period[0] = t2-t1;
		info->box.period[1] = t4-t3;
		info->box.period[2] = t6-t5;
		//Atoms		
		info->type.clear();
		std::getline(is, temp);
		
		coordlist_t x_temp(n);
		std::vector<std::string> type_temp(n);
		
		for (int i=0; i<n; i++) {
			is >> particle_id >> tempi >> t1 >> t2 >> t3;
			std::ostringstream tmp;
			tmp << tempi;
			type_temp[particle_id-1]  = tmp.str();
			coord_t xi;
			xi.push_back(t1*info->box.period[0] + info->box.boxlo[0]);
			xi.push_back(t2*info->box.period[1] + info->box.boxlo[1]);
			xi.push_back(t3*info->box.period[2] + info->box.boxlo[2]);
			
			//make sure the particle is within the box			
			for(int k=0; k<3; k++)
			{
				if(xi[k] > info->box.boxhi[k])
					xi[k] -= info->box.period[k];
				if(xi[k] < info->box.boxlo[k])
					xi[k] += info->box.period[k];
				
			}
			x_temp[particle_id-1] = xi;
		}
		for (int i=0; i<n; i++)
		{
			x.push_back(x_temp[i]);
			info->type.push_back(type_temp[i]);
		}
		is >> spc;
	}
	
	
	std::ostream& operator << (std::ostream& os, box_t& box)
	{
		os << box.boxlo << "\n" << box.boxhi << "\n" << box.period << "\n";
		return os;
	}

	std::ostream& operator << (std::ostream& os, coordlist_t& x)
	{
		for (unsigned int i=0; i<x.size(); i++) {
			os << x[i] <<"\n";
		}
		return os;
	}
	
	std::ostream& operator << (std::ostream& os, coord_t& x)
	{
		unsigned int sm1 = x.size()-1;
		for (unsigned int k=0; k<x.size(); k++) {
			os << x[k];
			if (k != sm1) {
				os <<"\t";
			}
		}
		return os;
	}
	
	std::ostream& operator << (std::ostream& os, shpdesc_t& sd)
	{
		unsigned int sm1 = sd.size()-1;
		for (unsigned int k=0; k<sd.size(); k++) {
			os << sd[k];
			if (k != sm1) {
				os <<"\t";
			}
		}
		return os;		
	}

	std::ostream& operator << (std::ostream& os, function1d_t& f)
	{
		for (function1d_t::iterator i=f.begin(); i!=f.end(); ++i) {
			os << i->first << "\t" << i->second << "\n";
		}
		return os;
	}


	std::istream& operator >> (std::istream& is, coordlist_t& x)
	{
		return is;
	}
	
	std::istream& operator >> (std::istream& is, coord_t& x)
	{
		return is;
	}
	
	std::istream& operator >> (std::istream& is, shpdesc_t& sd)
	{
		return is;
	}
}
