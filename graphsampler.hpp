#ifndef GRAPHSAMPLER_HPP_
#define GRAPHSAMPLER_HPP_

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "myutil.hpp"


namespace ublas = boost::numeric::ublas;
typedef ublas::symmetric_matrix<double> smat;
typedef ublas::matrix<double> mat;
typedef ublas::vector<double> vec;

using namespace std;

class Structure;
class Data;

class Data
{
	friend class Structure;
	friend bool readData(const char* filename, ublas::matrix<double>& data, vector< vector<MySizeType> >& missingIndex);
	
	public:
		double gamma;
		MySizeType getNumNodes() { return n; };
		Data(const char* );
		double getLocalScore(MySizeType, Structure&);
		void imputeData(MySizeType, Structure&);
		double test;

	private:
		MySizeType n;
		double m;
		double alpha;
		smat T;		
		double nu;
		double lcwishart(MySizeType, MySizeType);
		double base_score;
		vector< vector<double> > lcwishart_cache;
		double arg2[4];
		mat data;
		vector< vector<MySizeType> > missingIndex;
		mat decompMat;
		smat totalS;
		mat totalM;
		mat T0;
};

class Structure
{
	friend class Data;
	
	public:
		double mTemperature;
		vector< vector<MySizeType> > mParents;
		bool moveMCMC(Data&);
		double getScore();
		vector<MySizeType> getMB(MySizeType);
		Structure() {};
		Structure(Data& su)
		:mTemperature(1), mParents(su.n), mNumNodes(su.n), mLocalScores(su.n), isVisited(su.n), mChildren(su.n) { 
			randomize(mNumNodes/2, su);
		};

	private:
		MySizeType mNumNodes;
		vector<double> mLocalScores;
		void randomize(MySizeType, Data&);
		bool isDescendant(MySizeType, MySizeType);
		bool isDescendSub(MySizeType, MySizeType);
		vector<bool> isVisited;
		vector< vector<MySizeType> > mChildren;
};
		

#endif /*GRAPHSAMPLER_HPP_*/
