/*
 *  HistogramConst.h
 *  glotzmath
 *
 *  Created by Aaron Keys on 1/26/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GLOTZILLA_CONSTHISTOGRAM_H
#define _GLOTZILLA_CONSTHISTOGRAM_H

#include <vector>
#include <iostream>

namespace smac
{
	/**
	@brief a histogram with a constant size
	@author Aaron Keys
	*/
	
	template <class _TYPE>
	class HistogramConst
	{
		public:			
			/**
			@brief Constructor taking two arguments
			@param min is the minimum of the histogram
			@param max is the max of the histogram
			*/
			HistogramConst(_TYPE min, _TYPE step, _TYPE max) :
				_min(min),
				_inc(step),
				_max(max),
				_show_warnings(false)
			{
				_inv_range = 1.0 / (max - min);
				_nbins = (int)((max-min)/_inc);
				_weight_of_bin = new double[_nbins];
				clear();
			}
			
			~HistogramConst()
			{
				delete _weight_of_bin;
			}

			/**
			@brief Copy constructor
			@param rhs is another histogram of equal size
			*/
			HistogramConst(HistogramConst<_TYPE> &rhs)
			{
			}
			
			/**
			@brief insert a datapoint into the histogram
			@param value is the value of the datapoint
			@sa insert
			*/
			void operator += (_TYPE value)
			{
				insert(value);
			}

			/**
			@brief removes a datapoint into the histogram
			@param value is the value of the datapoint
			@sa remove
			*/
			void operator -= (_TYPE value)
			{
				remove(value);
			}
			
			/**
			@brief insert a datapoint into the histogram
			@param value is the value of the datapoint
			*/
			void insert(_TYPE value)
			{
				insert(value, 1.0);
			}
			
			/**
			@brief insert a datapoint into the histogram
			@param value is the value of the datapoint
			@param weight is the weight of the datapoint
			*/
			void insert(_TYPE value, double weight)
			{
				int bin = (int)((value * _inv_range) * _nbins);
				if(bin >=0 && bin < _nbins)
					_weight_of_bin[bin] += weight;
				else if(_show_warnings)
					std::cerr << "Warning: histogram value out of range\n";
			}


			/**
			@brief removes a datapoint into the histogram
			@param value is the value of the datapoint
			*/
			void remove(_TYPE value)
			{
				remove(value, 1.0);
			}
			
			/**
			@brief removes a datapoint into the histogram
			@param value is the value of the datapoint
			*/
			void remove(_TYPE value, double weight)
			{
				int bin = (int)((value * _inv_range) * _nbins);
				if(bin >=0 && bin < _nbins)
					_weight_of_bin[bin] -= weight;
				else if(_show_warnings)
					std::cerr << "Warning: histogram value out of range\n";
			}

			/**
			@brief indexing operator
			@param i is the index to be returned
			*/
			double & operator [](int i) const
			{
				return _weight_of_bin[i];
			}
			
			/**
			@brief clears the contents of the bins
			@sa Reset()
			*/
			void clear()
			{
				for(int i=0; i<_nbins; i++)
					_weight_of_bin[i] = 0.0;
			}
			
			/**
			@brief clears the contents of the bins
			@sa clear()
			*/
			void reset()
			{
				clear();
			}

			/**
			@brief prints the histogram to the stdout
			*/
			void print()
			{
				print(std::cout);
			}
			
			/**
			@brief 
			*/
			void print(std::ostream &os)
			{
				for(int i=0; i<_nbins; i++)
					os << (_TYPE)(_min + (_inc*i + (_inc*(i+1)))*0.5) <<"\t"<< _weight_of_bin[i] << std::endl; 
			}
			
			/**
			@brief Returns the contents of the bin at index i as a std::pair
			*/
			std::pair<_TYPE, _TYPE> getBin(unsigned int i)
			{
				if (i>= (unsigned int)_nbins) {
					std::cerr << "Error in glotzmath::HistogramDynamic: " <<
								 "Attempt to access out of bounds index";
				}
				
				std::pair<_TYPE, _TYPE> x;
				x.first = (_TYPE)(_min + ((_inc*i)+(_inc*(i+1)))*0.5);
				x.second = _weight_of_bin[i];
				return x;
			}

			/**
			@brief 
			*/
			int getNumberOfBins() const
			{
				return _nbins;
			}

		protected:
			_TYPE _min;															///< Minimum value
			_TYPE _inc;															///< Increment
			_TYPE _max;															///< Maximum value
			_TYPE _inv_range;													///< 1 over the range of values
			int _nbins;															///< Number of bins
			bool _show_warnings;												///< Indicates whether warnings should be displayed
			double *_weight_of_bin;												///< Bin weight (number of datapoints in each bin)
	};
}

#endif
