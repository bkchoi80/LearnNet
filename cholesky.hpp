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

#ifndef CHOLESKY_HPP_
#define CHOLESKY_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/operation.hpp>


namespace boost { namespace numeric { namespace ublas {

  // Cholesky factorization
  template<class M>
  void cholesky_factorize (M &m) {
    typedef M matrix_type;
    typedef typename M::size_type size_type;
    typedef typename M::value_type value_type;

    BOOST_UBLAS_CHECK (m.size1() == m.size2(), external_logic("Cholesky decomposition is only valid for a square, positive definite matrix."));

    size_type size = m.size1();
    vector<value_type> d(size);

    for (size_type i = 0; i < size; ++ i) {
      matrix_row<M> mri (row (m, i));
      for (size_type j = i; j < size; ++ j) {
        matrix_row<M> mrj (row (m, j));

        value_type elem = m(j,i) - inner_prod(project(mri,range(0,i)), project(mrj,range(0,i)));

        if (i == j) {
        	BOOST_UBLAS_CHECK (elem > 0.0, external_logic("The matrix is not positive definite."));
            d(i) = sqrt(elem);
        }
        else {
          m(j,i) = elem / d(i);
        }
      }
    }
    for (size_type i = 0; i < size; ++ i) {
      m(i,i) = d(i);
    }
  }
  
  // Cholesky update. Insert into the second to the bottom row.
  template<class M>
  void cholesky_update (M &m, typename M::size_type size, bool atBottom = false) {
    typedef typename M::size_type size_type;
	
    BOOST_UBLAS_CHECK (m.size1() == m.size2() && m.size1() > size, external_logic("Wrong dimesions."));

    for(size_type i=0 ; i<size ; i++)
    	m(size, i) = (m(size, i) - inner_prod(project(row(m, i), range(0, i)), project(row(m, size), range(0, i)))) / m(i, i);
    m(size, size) = m(size, size) - inner_prod(project(row(m, size), range(0, size)), project(row(m, size), range(0, size)));    
    BOOST_UBLAS_CHECK (m(size, size) > 0.0, external_logic("The matrix is not positive definite."));
    m(size, size) = sqrt(m(size, size));

    if(!atBottom)
    	cholesky_swap(m, size + 1, size - 1);
  }
  
  template<class M>
  void cholesky_downdate (M &m, typename M::size_type size, typename M::size_type dr) {
      typedef typename M::size_type size_type;
      typedef typename M::value_type value_type;

      BOOST_UBLAS_CHECK (m.size1() == m.size2() && m.size1() >= size && size > dr, external_logic("Dimension error!"));

      for(size_type i = dr; i < size - 1; i++)
    	  cholesky_swap(m, size, i);
  }
  
  template<class M>
  void cholesky_swap (M &m, typename M::size_type size, typename M::size_type sr) {
      typedef typename M::size_type size_type;
      typedef typename M::value_type value_type;
      
      BOOST_UBLAS_CHECK (m.size1() == m.size2() && m.size1() >= size && size > sr + 1, external_logic("Dimension error!"));

      size_type idx;
      
      value_type oldRow1 = m(sr, sr);
      value_type oldRow2[2] = {m(sr + 1, sr), m(sr + 1, sr + 1)};
      value_type newRow1 = sqrt(oldRow2[0]*oldRow2[0] + oldRow2[1]*oldRow2[1]);
      value_type newRow2[2];
      newRow2[0] = oldRow1 * oldRow2[0] / newRow1;
      newRow2[1] = sqrt(oldRow1*oldRow1 - newRow2[0]*newRow2[0]);
      
      for(idx = sr + 2; idx < size ; idx++){
    	  value_type newElem = (m(idx, sr)*oldRow2[0] + m(idx, sr + 1)*oldRow2[1]) / newRow1;
    	  m(idx, sr + 1) = (m(idx, sr)*oldRow1 - newElem*newRow2[0]) / newRow2[1];
    	  m(idx, sr) = newElem;
      }
      
      for(idx = 0; idx < sr; idx++)
    	  std::swap(m(sr, idx), m(sr + 1, idx));

      m(sr, sr) = newRow1;
      m(sr + 1, sr) = newRow2[0];
      m(sr + 1, sr + 1) = newRow2[1];
  }

}}} 


#endif /*CHOLESKY_HPP_*/
