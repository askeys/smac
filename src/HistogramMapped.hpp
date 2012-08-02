/*
 *  HistogramMapped.h
 *  glotzilla
 *
 *  Created by askeys on 7/17/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GLOTZILLA_MAPHISTOGRAM_H
#define _GLOTZILLA_MAPHISTOGRAM_H

#include <map>
#include <vector>
#include <iostream>

namespace smac
{
	/**
	@brief A histogram 
	*/
	template <class _TYPE>
	class HistogramMapped
	{
		private:
		
			typedef std::map <_TYPE, double> MapT; 
			MapT _key_vs_count;

		public:
			
			HistogramMapped()
			{
			}
			
			void clear()
			{
				_key_vs_count.clear();
			}
						
			HistogramMapped & operator = (const HistogramMapped &rhs)
			{
				_key_vs_count = rhs;
				return *this;
			}
			
			
			void insert(_TYPE value, double weight=1)
			{
				if(_key_vs_count.find(value) == _key_vs_count.end()) {
					_key_vs_count[value] = 0;
				}
				_key_vs_count[value] += weight;
			}
			
			void remove(_TYPE value, double weight=1)
			{
				if(_key_vs_count.find(value) == _key_vs_count.end()) {
					_key_vs_count.erase(value);
				}				
			}
			
			std::vector<_TYPE> getX()
			{
				std::vector<_TYPE> x;
				typename std::map< _TYPE, double >::iterator xf;
				
				for (xf = _key_vs_count.begin(); xf!= _key_vs_count.end(); ++xf) {
					x.push_back(xf->first);
				}
				
				return x;
			}
			
			std::vector<double> getF()
			{
				std::vector<double> f;
				typename std::map< _TYPE, double >::iterator xf;
				
				for (xf = _key_vs_count.begin(); xf!= _key_vs_count.end(); ++xf) {
					f.push_back(xf->second);
				}
				
				return f;
			}
			
			std::vector< std::pair<_TYPE,double> > getFunction()
			{
				std::vector< std::pair<_TYPE,double> > fx;
				typename std::map< _TYPE, double >::iterator xf;
				
				for (xf = _key_vs_count.begin(); xf!= _key_vs_count.end(); ++xf) {
					fx.push_back(*xf);
				}
				return fx;		
			}					
			
			void print(std::ostream &os)
			{
				typename MapT::iterator i = _key_vs_count.begin();
				typename MapT::iterator e = _key_vs_count.end();
				
				for( i; i != e; ++i ) {
					os << i->first <<"\t"<< i->second <<"\n";
				}
			}

			int getNumberOfBins()
			{
				return _key_vs_count.size();
			}
	};
}

#endif
