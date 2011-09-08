/*
 *	Copyright (c) 2009 Bokyung Choi <bkchoi@stanford.edu>
 *
 *	This file is part of LearnNet.
 *
 *	LearnNet is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	LearnNet is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with Foobar.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef MYUTIL_HPP_
#define MYUTIL_HPP_

#include <vector>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>


namespace myutil{
	void seed(size_t);
	size_t rand(size_t);
	double rand();
	bool readDelim(const char* const, std::vector< std::vector<double> >&);
	
	
	template <class T>
	void printVector(const std::vector<T>& vec)
	{
		for(typename std::vector<T>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
			std::cout << *iter << "\t";
		std::cout << std::endl;
	}
	
	template <class T>
	class Ranker
	{
	private:
		const std::vector<T>& vec;

	public:
		Ranker(const std::vector<T>& vec_) : vec(vec_) { }
		bool operator()(size_t idx1, size_t idx2) const
		{ 
			return(vec[idx1] < vec[idx2]);
		}
		void getOrder(std::vector<size_t>& order)
		{
			order.resize(vec.size());
			for(size_t idx = 0; idx < order.size(); idx++)
				order[idx] = idx;
			std::sort(order.begin(), order.end(), *this);
		}
	};	
	
	template <class T>
	bool sample(std::vector<T>& vec, size_t size)
	{
		if(size > vec.size())
			return false;
		
		if(size <= vec.size() / 2) {
			for(size_t sampleNum = 0; sampleNum < size; sampleNum++)
				std::iter_swap(vec.begin() + sampleNum, 
						vec.begin() + rand(vec.size() - sampleNum) + sampleNum);
		}
		if(size > vec.size() / 2) {
			for(size_t discardNum = 0; discardNum < vec.size() - size; discardNum++)
				std::iter_swap(vec.end() - 1 - discardNum, 
						vec.begin() + rand(vec.size() - discardNum));
		}
		vec.resize(size);

		return true;
	}
	
	template <class T>
	void getRandomPermutation(std::vector<T>& vec)
	{
		for(size_t elem = 0; elem < vec.size(); elem++)
			std::iter_swap(vec.begin() + elem, 
					vec.begin() + rand(vec.size() - elem) + elem);
	}
}


#endif /*MYUTIL_HPP_*/
