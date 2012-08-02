/*
 *  HistogramDynamic.h
 *  Data
 *
 *  Created by Aaron Keys on 1/26/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GLOTZILLA_DYNAMIC_HISTOGRAM_H
#define _GLOTZILLA_DYNAMIC_HISTOGRAM_H
#include <deque>
#include <iostream>

#include <cstdlib>

namespace smac
{
	/**
	@brief A dynamically sized histogram
	@author Aaron Keys
	@ingroup math
	@code
	#include <cstdlib>															//gives us the drand48() random number generator used in the example
	#include <iostream>															//gives us std::cerr used in the example

	int main()
	{
		glotzilla::HistogramDynamic <double> h(0.1);							//histogram "h" binning values of type double with binsize 0.1
		h.insert(drand48());													//insert a random double [0,1]
		h += drand48();															//another way to add a value to h
		h.insert(5);															//doesn't matter that it is out of bounds, h will resize to accomadate it
		h.Remove(drand48());													//de-increment the bin containing a number
		h -= drand48();															//another way to do it
		h.print();																//print to stdout
		h.print(std::cerr);														//print to stderr
	}

	@endcode
	*/

	template <class T>
	class HistogramDynamic
	{
		public:
			/**
			@brief Default constructor
			*/
			HistogramDynamic() :
				_inc(0.01),
				_is_firsttime(true)
			{
			}

			/**
			@brief Default constructor taking bin increment
			*/
			HistogramDynamic(T inc) :
				_inc(inc),
				_is_firsttime(true)
			{
			}

			/**
			@brief Copy constructor
			*/
			HistogramDynamic(const HistogramDynamic &rhs)
			{
				_weight_of_bin = rhs._weight_of_bin;
				_min = rhs._min;
				_max = rhs._max;
				_inc = rhs._inc;
				_is_firsttime = rhs._is_firsttime;
			}

			/**
			@brief Returns the size of each bin
			*/
			T getBinSize() const
			{
				return _inc;
			}

			/**
			@brief Addition function
			*/
			void combineWith(const HistogramDynamic& rhs)
			{
				if (_inc != rhs._inc) {
					std::cerr << "Error: glotzmath::HistogramDynamic:\n"<<
								"Attempt to combine histograms with different bin size\n";
					exit(-1);
				}

				HistogramDynamic temp(rhs);
				if (_max > temp._max) {
					temp.insert(_max-_inc/2.0, 0);
				}
				if (temp._max > _max) {
					insert(temp._max-temp._inc/2.0, 0);
				}
				if (_min < temp._max) {
					temp.insert(_min+_inc/2.0, 0);
				}
				if (temp._min < _min) {
					insert(temp._min+temp._inc/2.0, 0);
				}

				for (unsigned int i=0; i<temp._weight_of_bin.size(); i++) {
					_weight_of_bin[i] += temp._weight_of_bin[i];
				}
			}

			/**
			@brief Sets the bin size
			*/
			void setBinSize(T size)
			{
				_inc = size;
				clear();
			}

			/**
			@brief Indexing operator
			*/
			double & operator [](int i)
			{
				return _weight_of_bin[i];
			}

			/**
			@brief inserts a datapoint into the histogram
			@value is the value to be inserted
			*/
			void insert(T value)
			{
				insert(value, 1.0);
			}

			/**
			@brief inserts a datapoint into the histogram
			@value is the value to be inserted
			@weight is the weight of the datapoint
			*/
			void insert(T value, double weight)
			{
				if (_is_firsttime) {
					createFirstBin(value, weight);
					return;
				}

				int bin = computeBinContaining(value);

				if (bin == TOO_BIG) {
					expandRight(value, weight);
				}
				else if (bin == TOO_SMALL) {
					expandLeft(value, weight);
				}
				else {
					_weight_of_bin[bin] += weight;
				}
			}

			/**
			@brief clears the contents of the histogram
			@sa Reset()
			*/
			void clear()
			{
				_weight_of_bin.clear();
				_is_firsttime = true;
			}

			/**
			@brief clears the contents of the histogram
			@sa clear()
			*/
			void reset()
			{
				clear();
			}

			/**
			@brief prints the histogram to stdout
			*/
			void print()
			{
				print(std::cout);
			}

			/**
			@brief prints the histogram to a stream
			@os is the stream to print to
			*/
			void print(std::ostream &os)
			{
				for (unsigned int i=0; i<_weight_of_bin.size(); i++) {
					os	<< (T)(_min + ((_inc*i)+(_inc*(i+1)))*0.5)
						<<"\t"<< _weight_of_bin[i] << std::endl;
				}
			}

			/**
			@brief Returns the number of bins in the histogram
			*/
			unsigned int getNumberOfBins() const
			{
				return _weight_of_bin.size();
			}

			/**
			@brief Returns the contents of the bin at index i as a std::pair
			*/
			std::pair<T, T> getBin(unsigned int i)
			{
				if (i>= _weight_of_bin.size()) {
					std::cerr << "Error in glotzmath::HistogramDynamic: " <<
								 "Attempt to access out of bounds index";
				}

				std::pair<T, T> x;
				x.first = (T)(_min + ((_inc*i)+(_inc*(i+1)))*0.5);
				x.second = _weight_of_bin[i];
				return x;
			}

			std::pair<T, T> max()
			{
				std::pair<T, T> max;
				max.second = 1e-32;

				for (unsigned int i=0; i<_weight_of_bin.size(); i++) {
					if (_weight_of_bin[i] > max.second) {
						max.first = (T)(_min + ((_inc*i)+(_inc*(i+1)))*0.5);
						max.second = _weight_of_bin[i];
					}
				}
				return max;
			}

		protected:

			/**
			@brief expands the histogram in the positive x direction
			@value is the value of the out-of-bounds data point
			@weight is the weight of the point
			*/
			void expandRight(T value, double weight)
			{
				for (T i = _max; i < value; i+=_inc)
				{
					_max = i + _inc;
					_weight_of_bin.push_back(0);
				}
				_weight_of_bin[_weight_of_bin.size()-1] += weight;
			}

			/**
			@brief expands the histogram in the negative x direction
			@value is the value of the out-of-bounds datapoint
			@weight is the weight of the point
			*/
			void expandLeft(T value, double weight)
			{
				for (T i = _min; i > value; i-=_inc) {
					_min = i - _inc;
					_weight_of_bin.push_front(0);
				}
				_weight_of_bin[0] += weight;
			}

			/**
			@brief computes the bin number containing a given datapoint
			@value is the value of the datapoint
			*/
			int computeBinContaining(T value)
			{
				if (value < _min) {
					return TOO_SMALL;
				}
				else if (value >= _max) {
					return TOO_BIG;
				}
				else {
					return (int)( (value - _min) / (_max-_min) *
						_weight_of_bin.size() );
				}
			}

			/**
			@brief Creates the initial bin in the histogram
			*/
			void createFirstBin(T value, double weight)
			{
				_weight_of_bin.push_back(weight);
				_min = value - _inc * .5;
				_max = value + _inc * .5;
				_is_firsttime = false;
			}

			std::deque <double> _weight_of_bin;		///< Dynamic sized container holding the weight associated with each bin
			T _min;									///< Minumum bound of the histogram
			T _max;									///< Maximum bound of the histogram
			T _inc;									///< Increment of the histogram
			enum{TOO_SMALL = -2, TOO_BIG = -1};		///< Holds values for indicating out-of-bounds indexing
			bool _is_firsttime;						///< Indicates whether the first bin exists

	};

	template <class TYPE>
	std::ostream& operator << (std::ostream& os, HistogramDynamic<TYPE>& h)
	{
		h.print(os);
		return(os);
	}

}

#endif
