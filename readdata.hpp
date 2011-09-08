#ifndef READDATA_HPP_
#define READDATA_HPP_

#include <boost/spirit.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "myutil.hpp"


using namespace boost::spirit;
namespace ublas = boost::numeric::ublas;

bool readData(const char*, ublas::matrix<double>&, vector< vector<MySizeType> >&);


#endif /*READDATA_HPP_*/
