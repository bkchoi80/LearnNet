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

#ifndef LEARNNET_HPP_
#define LEARNNET_HPP_

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include "myutil.hpp"


enum ScoreType {MLSCORE, BDEUSCORE};

class Dag
{
public:
	Dag() : mNumNodes(0), mNumEdges(0) {}
	Dag(size_t numNodes) {setNumNodes(numNodes);}
	size_t getNumNodes() const {return mNumNodes;}
	size_t getNumEdges() const {return mNumEdges;}
	const std::vector<size_t>& getParents(size_t node) const {
		assert(node < mNumNodes);
		return mParents[node];
	}
	const std::vector<size_t>& getOrder() {
		updateOrder();
		return mOrder;
	}
	bool addParent(size_t parent, size_t node, bool cycleCheck = true);
	bool removeParent(size_t parent, size_t node);
	void clearParents(size_t node) {mNumEdges -= mParents[node].size(); mParents[node].clear();}
	bool isDescendant(size_t node1, size_t node2) const;
	inline int findParent(size_t parent, size_t node) const;
	void findVStruct(size_t node, std::vector<size_t>& result) const;
	void randomize();
	void setNumNodes(size_t);
	void updateOrder();

protected:
	size_t mNumNodes;
	size_t mNumEdges;
	std::vector< std::vector<size_t> > mParents;
	std::vector<size_t> mOrder;
	bool mOrderUpdated;

private:
	mutable std::vector<bool> mVisited;
	bool isDescendantSub(size_t, size_t) const;
};

std::ostream& operator <<(std::ostream &os, const Dag &dag);

class Data
{
public:
	Data(const std::vector< std::vector<double> >& covariance, double dataSize);
	boost::numeric::ublas::symmetric_matrix<double> mCovariance;
	double mDataSize;
	size_t mNumNodes;
	double mNu;
	double mAlpha;
	double mBaseScore[2];
	std::vector<double> mBDeuCoeffs;
};

class OrderSearch : public Dag
{
public:
	OrderSearch() : Dag(), mpData(0) {}
	OrderSearch(const Data* pData, double gamma = 0, ScoreType scoreType = MLSCORE) {
		attachData(pData, gamma, scoreType);
	}
	void attachData(const Data* pData, double gamma = 0, ScoreType scoreType = MLSCORE);
	double getScore() const {return mTotalScore;}
	void randomize();
	const std::vector<double>& getLocalScores() const {
		return mLocalScores;
	}
	double getGamma() const {
		return mGamma;
	}
	double getNextGamma() const {
		return mNextGamma;
	}
	const Dag doTabuSearch(size_t numIter, size_t tabuPeriod, size_t restartPeriod, bool verbose = false, size_t numBestGraphs = 10, double epsilonScore = 0.00001);
	void moveToNextGamma() {
		if(mNextGamma > 0)
			setGamma(mNextGamma);
	}


private:
	double mGamma;
	const Data *mpData;
	ScoreType mScoreType;
	bool mOrderEvaluated;
	double mNextGamma;
	double mTotalScore;
	std::vector<double> mLocalScores;
	mutable std::vector< boost::numeric::ublas::matrix<double> > mCholCache;
	mutable std::vector< boost::numeric::ublas::matrix<double> > mSwapOrderCache;
	void setGamma(double);
	void setScoreType(ScoreType scoreType);
	bool addParent(size_t parent, size_t node, bool cycleCheck = true);
	bool removeParent(size_t parent, size_t node);
	inline void clearParents(size_t node);
	void swapOrder(size_t orderNum);
	double scoreSwap(size_t orderNum);
	void evaluateOrder();
	inline double computeScore(size_t) const;
	void findOptimalParents(size_t orderNum);
};

#endif /*LEARNNET_HPP_*/
