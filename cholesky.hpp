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
  bool cholesky_factorize (M &m) {
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
          if (elem <= 0.0) {
            // matrix after rounding errors is not positive definite
            return false;
          }
          else {
            d(i) = sqrt(elem);
          }
        }
        else {
          m(j,i) = elem / d(i);
        }
      }
    }
    // put the diagonal back in
    for (size_type i = 0; i < size; ++ i) {
      m(i,i) = d(i);
    }
    
    // decomposition succeeded
    return true;
  }
  
  // Cholesky factorization with submatrix
  template<class M>
  bool cholesky_factorize (M &m, unsigned size) {
    typedef M matrix_type;
    typedef typename M::size_type size_type;
    typedef typename M::value_type value_type;

    BOOST_UBLAS_CHECK ((m.size1() >= size) && (m.size2() >= size), external_logic("Cholesky decomposition is only valid for a square, positive definite matrix."));

    vector<value_type> d(size);

    for (size_type i = 0; i < size; i++) {
      matrix_row<M> mri (row (m, i));
      for (size_type j = i; j < size; j++) {
        matrix_row<M> mrj (row (m, j));

        value_type elem = m(j,i) - inner_prod(project(mri,range(0,i)), project(mrj,range(0,i)));

        if (i == j) {
          if (elem <= 0.0) {
            // matrix after rounding errors is not positive definite
            return false;
          }
          else {
            d(i) = sqrt(elem);
          }
        }
        else {
          m(j,i) = elem / d(i);
        }
      }
    }
    // put the diagonal back in
    for (size_type i = 0; i < size; i++) {
      m(i,i) = d(i);
    }
    
    // decomposition succeeded
    return true;
  }
  
  // Cholesky update. Insert into the second to the bottom row.
  template<class M>
  bool cholesky_update (M &m, typename M::size_type size, vector<typename M::value_type> &nr) {
    typedef typename M::size_type size_type;

    for(size_type i=0 ; i<size-1 ; i++)
    	nr(i) = ( nr(i) - inner_prod(project(row(m, i), range(0, i)), project(nr, range(0, i))) ) / m(i, i);
    
    // matrix after rounding errors is not positive definite
    if(nr(size-1) - pow(norm_2( project(nr, range(0, size-1)) ), 2) <= 0.0)
    	return(false);
    
    nr(size-1) = sqrt( nr(size-1) - pow(norm_2( project(nr, range(0, size-1)) ), 2) );
    nr(size) = ( nr(size) - inner_prod(project(nr, range(0, size-1)), project(row(m, size-1), range(0, size-1))) ) / nr(size-1);
    nr(size+1) = sqrt( pow(m(size-1, size-1), 2) - pow(nr(size), 2) );

    // update succeeded
    return true;
  }
  
  // Cholesky update. Insert into the second to the bottom row.
  template<class M>
  bool cholesky_downdate (M &m, typename M::size_type size, M &update_temp, typename M::size_type dr) {
    typedef typename M::size_type size_type;
    typedef typename M::value_type value_type;

    BOOST_UBLAS_CHECK (size > dr, external_logic("Dimension error!"));

    for(size_type i=dr ; i<size-1 ; i++){
    	matrix_row<M> mri(m, i+1);
    	matrix_row<M> nmri(update_temp, i);
    	for(size_type j=i ; j<size-1 ; j++){
 			matrix_row<M> mrj(m, j+1);
 			matrix_row<M> nmrj(update_temp, j);
 			value_type oldval = inner_prod(project(mri, range(dr, i+2)), project(mrj, range(dr, i+2)));
    		oldval -= inner_prod(project(nmri, range(dr, i)), project(nmrj, range(dr, i)));
    		if(i==j){
    			// matrix after rounding errors is not positive definite
    			if(oldval <= 0.0)
 	   				return(false);
    			
    			update_temp(j, i) = sqrt(oldval);
    		}
    		else
    			update_temp(j, i) = oldval / update_temp(i, i);
    	}
    }
    
    // downdate succeeded
    return true;
  }

}}} 


#endif /*CHOLESKY_HPP_*/
