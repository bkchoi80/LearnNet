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

#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

#include "myutil.hpp"
#include "cholesky.hpp"
#include "learnnet.hpp"


using namespace std;

ostream& operator <<(ostream &os,const Dag &dag)
{
	for(size_t nodeNum = 0; nodeNum < dag.getNumNodes(); nodeNum++) {
		os << "Node " << nodeNum << ":  ";
		for(size_t parentNum = 0; parentNum < dag.getParents(nodeNum).size(); parentNum++)
			os << "\t" << dag.getParents(nodeNum)[parentNum];
		os << endl;
	}
	os << endl;

	return os;
}

void Dag::setNumNodes(size_t numNodes)
{
	mParents.resize(numNodes);
	for(size_t nodeNum = 0; nodeNum < mParents.size(); nodeNum++)
		mParents[nodeNum].clear();
	mVisited.resize(numNodes);
	mOrder.resize(numNodes);
	mOrderUpdated = false;
	mNumNodes = numNodes;
	mNumEdges = 0;
}

inline int Dag::findParent(size_t parent, size_t node) const
{
	vector<size_t>::const_iterator iter = find(mParents[node].begin(), mParents[node].end(), parent);
	return (iter == mParents[node].end()) ? -1 : (int)(iter - mParents[node].begin());
}

void Dag::findVStruct(size_t node, vector<size_t>& result) const
{
	result.clear();
	result.resize(mParents[node].size(), 1);

	for(size_t parentNum1 = 0; parentNum1 < mParents[node].size(); parentNum1++) {
		for(size_t parentNum2 = parentNum1 + 1; parentNum2 < mParents[node].size(); parentNum2++) {
			if(result[parentNum1] || result[parentNum2] == 0)
				continue;
			if(findParent(mParents[node][parentNum1], mParents[node][parentNum2]) < 0)
				result[parentNum1] = result[parentNum2] = 0;
			else if(findParent(mParents[node][parentNum2], mParents[node][parentNum1]) < 0)
				result[parentNum1] = result[parentNum2] = 0;
		}
	}

	for(size_t parentNum = 0; parentNum < mParents[node].size(); parentNum++)
		if(result[parentNum])
			result.push_back(mParents[node][parentNum]);

	result.erase(result.begin(), result.begin() + mParents[node].size());
}

bool Dag::addParent(size_t parent, size_t node, bool cycleCheck)
{
	assert(node < mNumNodes && parent < mNumNodes);

	if(cycleCheck) {
		if(findParent(parent, node) != -1 || isDescendant(parent, node))
			return false;
	}
	else
		if(findParent(parent, node) != -1)
			return false;

	mParents[node].push_back(parent);
	mNumEdges++;
	mOrderUpdated = false;
	return true;
}

bool Dag::removeParent(size_t parent, size_t node)
{
	assert(node < mNumNodes && parent < mNumNodes);

	int found = findParent(parent, node);
	if(found != -1) {
		mParents[node].erase(mParents[node].begin() + found);
		mNumEdges--;
		return true;
	}
	else
		return false;
}

void Dag::randomize()
{
	for(size_t orderNum = 0; orderNum < mOrder.size(); orderNum++)
		mOrder[orderNum] = orderNum;
	myutil::getRandomPermutation(mOrder);

	for(size_t parentNum = 0; parentNum < mParents.size(); parentNum++)
		mParents[parentNum].clear();
	mNumEdges = 0;

	size_t node1, node2;
	while(mNumEdges < mNumNodes / (myutil::rand(2) + 1)) {
		node1 = myutil::rand(mNumNodes);
		node2 = myutil::rand(mNumNodes-1);
		if(node2 >= node1)
			node2++;
		if(node2 < node1)
			std::swap(node1, node2);

		addParent(mOrder[node1], mOrder[node2], true);
	}

	mOrderUpdated = true;
}

bool Dag::isDescendant(size_t node1, size_t node2) const
{
	for(size_t nodeNum = 0; nodeNum < mNumNodes; nodeNum++)
		mVisited[nodeNum] = 0;

	return isDescendantSub(node1, node2);
}

bool Dag::isDescendantSub(size_t node1, size_t node2) const
{
	size_t i;

	if(mVisited[node1])
		return(0);
	else
		mVisited[node1] = 1;

	for(i = 0; i < mParents[node1].size(); i++)
		if(mParents[node1][i] == node2)
			return(1);

	for(i=0 ; i<mParents[node1].size() ; i++)
		if(isDescendantSub(mParents[node1][i], node2))
			return(1);

	return(0);
}

void Dag::updateOrder()
{
	if(mOrderUpdated)
		return;

	size_t numEdges = mNumEdges;
	vector< vector<size_t> > parents(mParents);
	vector<bool> inserted(mNumNodes, false);

	for(size_t orderNum = 0; orderNum < mOrder.size(); orderNum++) {
		for(size_t nodeNum1 = 0; nodeNum1 < mParents.size(); nodeNum1++) {
			if(!inserted[nodeNum1] && mParents[nodeNum1].size() == 0) {
				mOrder[orderNum] = nodeNum1;
				for(size_t nodeNum2 = 0; nodeNum2 < mParents.size(); nodeNum2++) {
					if(!inserted[nodeNum2])
						removeParent(nodeNum1, nodeNum2);
				}
				inserted[nodeNum1] = true;
				break;
			}
		}
	}

	mParents = parents;
	mNumEdges = numEdges;
	mOrderUpdated = true;
}

Data::Data(const std::vector< std::vector<double> >& covariance, double dataSize)
{
	assert(covariance.size() > 0 && covariance.size() == covariance[0].size());
	mCovariance.resize(covariance.size());
	for(size_t rowNum = 0; rowNum < covariance.size(); rowNum++)
		for(size_t colNum = 0; colNum <= rowNum; colNum++)
			mCovariance(rowNum, colNum) = covariance[rowNum][colNum];

	mDataSize = dataSize;
	mNumNodes = covariance.size();
	mNu = 1.0;
	mAlpha = mNumNodes;
	mBDeuCoeffs.resize(mNumNodes);

	for(size_t i = 0; i < mBDeuCoeffs.size(); i++)
		mBDeuCoeffs[i] = lgamma((mAlpha + mDataSize - i) / 2.0) - lgamma((mAlpha - i) / 2.0);

	mBaseScore[MLSCORE] = -mDataSize / 2.0 * log(2.0 * M_PI * exp(1.0));
	mBaseScore[BDEUSCORE] = -mDataSize / 2.0 * log(M_PI)
							+ log(mNu / (mNu + mDataSize)) / 2.0
							- (mAlpha + mDataSize) / 2.0 * log(mDataSize);
}

bool OrderSearch::addParent(size_t parent, size_t node, bool cycleCheck)
{
	using namespace boost::numeric::ublas;

	if(!Dag::addParent(parent, node, cycleCheck))
		return false;

	if(mCholCache[node].size1() < mParents[node].size() + 1)
		mCholCache[node].resize(mParents[node].size() + 1, mParents[node].size() + 1);
	for(size_t parentNum = 0; parentNum < mParents[node].size() - 1; parentNum++)
		mCholCache[node](mParents[node].size(), parentNum) = mpData->mCovariance(parent, mParents[node][parentNum]);
	mCholCache[node](mParents[node].size(), mParents[node].size() - 1) = mpData->mCovariance(parent, node);
	mCholCache[node](mParents[node].size(), mParents[node].size()) = mpData->mCovariance(parent, parent);

	cholesky_update(mCholCache[node], mParents[node].size());

	double oldScore = mTotalScore;
	mTotalScore -= mLocalScores[node];
	mLocalScores[node] = computeScore(node);
	mTotalScore += mLocalScores[node];

	if(mTotalScore > oldScore) {
		double localNextGamma = mGamma + (mTotalScore - oldScore) * 2 / log(mpData->mDataSize);
		mNextGamma = (mNextGamma > localNextGamma || mNextGamma < 0) ? localNextGamma : mNextGamma;
	}

	mOrderEvaluated = false;

	return true;
}

bool OrderSearch::removeParent(size_t parent, size_t node)
{
	using namespace boost::numeric::ublas;

	assert(node < mNumNodes && parent < mNumNodes);

	int found = findParent(parent, node);
	if(found == -1)
		return false;
	else {
		mParents[node].erase(mParents[node].begin() + found);
		mNumEdges--;
	}

	cholesky_downdate(mCholCache[node], mParents[node].size() + 2, found);

	mTotalScore -= mLocalScores[node];
	mLocalScores[node] = computeScore(node);
	mTotalScore += mLocalScores[node];

	mOrderEvaluated = false;

	return true;
}

void OrderSearch::attachData(const Data* pData, double gamma, ScoreType scoreType)
{
	using namespace boost::numeric::ublas;

	mpData = pData;
	Dag::setNumNodes(pData->mNumNodes);
	mLocalScores.resize(mNumNodes);
	mCholCache.resize(mNumNodes);
	mSwapOrderCache.resize(mNumNodes);
	mScoreType = scoreType;

	for(size_t nodeNum = 0; nodeNum < mSwapOrderCache.size(); nodeNum++)
		mSwapOrderCache[nodeNum].resize(mNumNodes, mNumNodes);

	mTotalScore = 0;
	mNextGamma = -1;
	for(size_t nodeNum = 0; nodeNum < mNumNodes; nodeNum++) {
		mCholCache[nodeNum].resize(1, 1);
		mCholCache[nodeNum](0, 0) = sqrt(mpData->mCovariance(nodeNum, nodeNum));
	}
	setGamma(gamma);

	mOrderEvaluated = false;
}

void OrderSearch::randomize()
{
	Dag::randomize();

	mTotalScore = 0;
	for(size_t nodeNum = 0; nodeNum < mNumNodes; nodeNum++) {
		mCholCache[nodeNum].resize(mParents[nodeNum].size() + 1, mParents[nodeNum].size() + 1);
		for(size_t rowNum = 0; rowNum < mParents[nodeNum].size(); rowNum++)
			for(size_t colNum = 0; colNum < mParents[nodeNum].size(); colNum++)
				mCholCache[nodeNum](rowNum, colNum) = mpData->mCovariance(mParents[nodeNum][rowNum], mParents[nodeNum][colNum]);
		for(size_t colNum = 0; colNum < mParents[nodeNum].size(); colNum++) {
			mCholCache[nodeNum](mParents[nodeNum].size(), colNum) = mpData->mCovariance(nodeNum, mParents[nodeNum][colNum]);
			mCholCache[nodeNum](colNum, mParents[nodeNum].size()) = mpData->mCovariance(mParents[nodeNum][colNum], nodeNum);
		}
		mCholCache[nodeNum](mParents[nodeNum].size(), mParents[nodeNum].size()) = mpData->mCovariance(nodeNum, nodeNum);

		cholesky_factorize(mCholCache[nodeNum]);
		mLocalScores[nodeNum] = computeScore(nodeNum);
		mTotalScore += mLocalScores[nodeNum];
	}
	mNextGamma = -1;

	mOrderEvaluated = false;
}

void OrderSearch::setGamma(double gamma)
{
	mGamma = gamma;
	mNextGamma = -1;
	mTotalScore = 0;
	for(size_t nodeNum = 0; nodeNum < mNumNodes; nodeNum++) {
		mLocalScores[nodeNum]  = computeScore(nodeNum);
		mTotalScore += mLocalScores[nodeNum];
	}

	mOrderEvaluated = false;
}

void OrderSearch::setScoreType(ScoreType scoreType)
{
	mScoreType = scoreType;
	mTotalScore = 0;
	mNextGamma = -1;
	for(size_t nodeNum = 0; nodeNum < mNumNodes; nodeNum++) {
		mLocalScores[nodeNum]  = computeScore(nodeNum);
		mTotalScore += mLocalScores[nodeNum];
	}

	mOrderEvaluated = false;
}

inline double OrderSearch::computeScore(size_t node) const
{
	double score = mpData->mBaseScore[mScoreType]
	            - mpData->mDataSize / 2.0 * log(mCholCache[node](mParents[node].size(), mParents[node].size()))
				- log(mpData->mDataSize) / 2.0 * mGamma * mParents[node].size();

	if(mScoreType == BDEUSCORE)
		score += log(mpData->mDataSize) / 2.0 * mParents[node].size()
				- mpData->mAlpha * log(mCholCache[node](mParents[node].size(), mParents[node].size()))
				+ mpData->mBDeuCoeffs[mParents[node].size()];

	return score;
}

void OrderSearch::clearParents(size_t node)
{
		Dag::clearParents(node);
		mCholCache[node](0, 0) = sqrt(mpData->mCovariance(node, node));

		mTotalScore -= mLocalScores[node];
		mLocalScores[node] = computeScore(node);
		mTotalScore += mLocalScores[node];

		mOrderEvaluated = false;
}

void OrderSearch::evaluateOrder()
{
	using namespace boost::numeric::ublas;

	if(mOrderEvaluated || mNumNodes == 1)
		return;

	updateOrder();

	for(size_t rowNum = 0; rowNum < mNumNodes; rowNum++)
		for(size_t colNum = 0; colNum < mNumNodes; colNum++)
			mSwapOrderCache[0](rowNum, colNum) = mpData->mCovariance(mOrder[rowNum], mOrder[colNum]);
	cholesky_factorize(mSwapOrderCache[0]);
	for(size_t matrixNum = 1; matrixNum < mSwapOrderCache.size(); matrixNum++)
		mSwapOrderCache[matrixNum] = mSwapOrderCache[0];
	for(size_t matrixNum = 0; matrixNum < mSwapOrderCache.size(); matrixNum++)
		cholesky_downdate(mSwapOrderCache[matrixNum], mSwapOrderCache[matrixNum].size1(), matrixNum);

	for(size_t orderNum = 0; orderNum < mOrder.size(); orderNum++)
		findOptimalParents(orderNum);

	mOrderEvaluated = true;
	mOrderUpdated = true;
}

void OrderSearch::findOptimalParents(size_t orderNum)
{
	if(orderNum == 0) {
		clearParents(mOrder[orderNum]);
		return;
	}

	vector<size_t> parentRank(orderNum);
	vector<double> residuals(orderNum);
	for(size_t parentNum = 0; parentNum < orderNum; parentNum++)
		residuals[parentNum] = -abs(mSwapOrderCache[parentNum](orderNum-1, orderNum-1));
	myutil::Ranker<double> ranker(residuals);
	ranker.getOrder(parentRank);

	size_t preserved = 0;
	while(preserved < min(orderNum, mParents[mOrder[orderNum]].size())) {
		if(mOrder[parentRank[preserved]] != mParents[mOrder[orderNum]][preserved])
			break;
		preserved++;
	}

	if(preserved < mParents[mOrder[orderNum]].size()) {
		for(size_t colNum = 0; colNum < preserved; colNum++)
			mCholCache[mOrder[orderNum]](preserved, colNum) = mCholCache[mOrder[orderNum]](mParents[mOrder[orderNum]].size(), colNum);

		mCholCache[mOrder[orderNum]](preserved, preserved) = pow(mCholCache[mOrder[orderNum]](mParents[mOrder[orderNum]].size(), mParents[mOrder[orderNum]].size()), 2);
		for(size_t colNum = preserved; colNum < mParents[mOrder[orderNum]].size(); colNum++)
			mCholCache[mOrder[orderNum]](preserved, preserved)
	     		= mCholCache[mOrder[orderNum]](preserved, preserved) + pow(mCholCache[mOrder[orderNum]](mParents[mOrder[orderNum]].size(), colNum), 2);
		mCholCache[mOrder[orderNum]](preserved, preserved) = sqrt(mCholCache[mOrder[orderNum]](preserved, preserved));
	}

	mNumEdges = mNumEdges - mParents[mOrder[orderNum]].size() + preserved;
	mParents[mOrder[orderNum]].erase(mParents[mOrder[orderNum]].begin() + preserved, mParents[mOrder[orderNum]].end());

	mTotalScore -= mLocalScores[mOrder[orderNum]];
	mLocalScores[mOrder[orderNum]] = computeScore(mOrder[orderNum]);
	mTotalScore += mLocalScores[mOrder[orderNum]];

	double bestScore;
	size_t unimproved = 0;
	for(size_t rankNum = mParents[mOrder[orderNum]].size(); rankNum < orderNum; rankNum++) {
		bestScore = mTotalScore;
		addParent(mOrder[parentRank[rankNum]], mOrder[orderNum], false);

		if(bestScore > mTotalScore)
			unimproved++;
		else {
			bestScore = mTotalScore;
			unimproved = 0;
		}

		if(unimproved > 2)
			break;
	}
	if(unimproved > 0)
		for(size_t removeNum = 0; removeNum < unimproved; removeNum++)
			removeParent(mParents[mOrder[orderNum]].back(), mOrder[orderNum]);

	if(mParents[mOrder[orderNum]].size() > 0) {
		size_t lastParent = mParents[mOrder[orderNum]].back();
		removeParent(lastParent, mOrder[orderNum]);
		addParent(lastParent, mOrder[orderNum]);
	}
}

void OrderSearch::swapOrder(size_t orderNum)
{
	using namespace boost::numeric::ublas;

	assert(orderNum >= 0 && orderNum < mNumNodes - 1);
	evaluateOrder();

	for(size_t matrixNum = 0; matrixNum < orderNum; matrixNum++)
		cholesky_swap(mSwapOrderCache[matrixNum], mNumNodes - 1, orderNum - 1);
	for(size_t matrixNum = orderNum + 2; matrixNum < mSwapOrderCache.size(); matrixNum++)
		cholesky_swap(mSwapOrderCache[matrixNum], mNumNodes - 1, orderNum);
	mSwapOrderCache[orderNum].swap(mSwapOrderCache[orderNum + 1]);

	swap(mOrder[orderNum], mOrder[orderNum+1]);

	if(Dag::findParent(mOrder[orderNum+1], mOrder[orderNum]) != -1)
		findOptimalParents(orderNum);
	findOptimalParents(orderNum+1);

	mOrderEvaluated = true;
	mOrderUpdated = true;
}

double OrderSearch::scoreSwap(size_t orderNum)
{
	using namespace boost::numeric::ublas;

	assert(orderNum >= 0 && orderNum < mNumNodes - 1);
	evaluateOrder();

	double score = mLocalScores[mOrder[orderNum]] + mLocalScores[mOrder[orderNum + 1]];

    for(size_t matrixNum = 0; matrixNum < orderNum; matrixNum++)
		cholesky_swap(mSwapOrderCache[matrixNum], orderNum + 1, orderNum - 1);
	mSwapOrderCache[orderNum].swap(mSwapOrderCache[orderNum + 1]);
	swap(mOrder[orderNum], mOrder[orderNum+1]);

	if(Dag::findParent(mOrder[orderNum+1], mOrder[orderNum]) != -1)
		findOptimalParents(orderNum);
	findOptimalParents(orderNum+1);

	score = mLocalScores[mOrder[orderNum]] + mLocalScores[mOrder[orderNum + 1]] - score;

	for(size_t matrixNum = 0; matrixNum < orderNum; matrixNum++)
		cholesky_swap(mSwapOrderCache[matrixNum], orderNum + 1, orderNum - 1);
	mSwapOrderCache[orderNum].swap(mSwapOrderCache[orderNum + 1]);
	swap(mOrder[orderNum], mOrder[orderNum+1]);

	if(Dag::findParent(mOrder[orderNum+1], mOrder[orderNum]) != -1)
		findOptimalParents(orderNum);
	findOptimalParents(orderNum+1);

	mOrderEvaluated = true;
	mOrderUpdated = true;

	return(getScore() + score);
}

const Dag OrderSearch::doTabuSearch(size_t numIter, size_t tabuPeriod, size_t restartPeriod, bool verbose, size_t numBestGraphs, double epsilonScore)
{
	size_t maxIdx;
	vector< vector<int> > tabuList(mNumNodes, vector<int> (mNumNodes, -1));
	vector<double> swapScores(mNumNodes - 1);
	vector<Dag> bestGraphs(numBestGraphs);
	vector<double> bestGraphNextGammas(numBestGraphs);
	vector<double> bestGraphScores(numBestGraphs);

	if(verbose)
		cout << endl << "Starting search with with penalty constant " << mGamma << "..." << endl;

	int restartCount = restartPeriod;
	double localMaxScore = mTotalScore;
	mNextGamma = -1;
	for(size_t iterNum = 0; iterNum < numIter; iterNum++) {
		if(restartCount < 0) {
			if(verbose) {
				cout << iterNum / 1000 << "k iterations done..." << endl;
				cout << "Local maximum: " << localMaxScore << "\t" << "Global Maximum: " << bestGraphScores.front() << endl << endl;
				cout << "Restarting..." << endl;
			}

			randomize();
			for(size_t orderNum = 0; orderNum < mNumNodes - 1; orderNum++)
				swapScores[orderNum] = scoreSwap(orderNum);

			for(size_t nodeNum1 = 0; nodeNum1 < mNumNodes; nodeNum1++)
				for(size_t nodeNum2 = 0; nodeNum2 < mNumNodes; nodeNum2++)
					tabuList[nodeNum1][nodeNum2] = -1;

			restartCount = restartPeriod;
			localMaxScore = mTotalScore;
		}

		maxIdx = swapScores.size();
		for(size_t scoreNum = 0; scoreNum < swapScores.size(); scoreNum++) {
			if(tabuList[mOrder[scoreNum]][mOrder[scoreNum+1]] < 0) {
				if(maxIdx == swapScores.size())
					maxIdx = scoreNum;
				else if(swapScores[scoreNum] > swapScores[maxIdx])
					maxIdx = scoreNum;
			}
		}

		swapScores[maxIdx] = mTotalScore;
		swapOrder(maxIdx);

		if(maxIdx > 0)
			swapScores[maxIdx - 1] = scoreSwap(maxIdx - 1);
		if(maxIdx < swapScores.size() - 1)
			swapScores[maxIdx + 1] = scoreSwap(maxIdx + 1);

		tabuList[mOrder[maxIdx]][mOrder[maxIdx+1]] = tabuPeriod;
		tabuList[mOrder[maxIdx+1]][mOrder[maxIdx]] = tabuPeriod;

		if(bestGraphs.back().getNumNodes() == 0 || mTotalScore > bestGraphScores.back() + epsilonScore) {
			for(size_t bestGraphNum = 0; bestGraphNum < numBestGraphs ; bestGraphNum++) {
				if(bestGraphs[bestGraphNum].getNumEdges() == 0) {
					bestGraphs[bestGraphNum] = (Dag) (*this);
					bestGraphScores[bestGraphNum] = mTotalScore;
					bestGraphNextGammas[bestGraphNum] = mNextGamma;

					break;
				}
				else if(mTotalScore > bestGraphScores[bestGraphNum] + epsilonScore) {
					for(size_t idx = bestGraphNum + 1; idx < numBestGraphs; idx++) {
						bestGraphs[idx] = bestGraphs[idx - 1];
						bestGraphScores[idx] = bestGraphScores[idx - 1];
						bestGraphNextGammas[idx] = bestGraphNextGammas[idx - 1];
					}
					bestGraphs[bestGraphNum] = (Dag) (*this);
					bestGraphScores[bestGraphNum] = mTotalScore;
					bestGraphNextGammas[bestGraphNum] = mNextGamma;

					break;
				}
				else if(mTotalScore > bestGraphScores[bestGraphNum] - epsilonScore)
					break;
			}
		}

		if(localMaxScore < mTotalScore)
			localMaxScore = mTotalScore;
		else
			restartCount--;

		for(size_t nodeNum1 = 0; nodeNum1 < mNumNodes; nodeNum1++)
			for(size_t nodeNum2 = 0; nodeNum2 < mNumNodes; nodeNum2++)
				if(tabuList[nodeNum1][nodeNum2] >= 0)
					tabuList[nodeNum1][nodeNum2]--;
	}


	size_t numEdges = bestGraphs.front().getNumEdges();
	double maxScore = bestGraphScores.front();
	double nextGamma = bestGraphNextGammas.front();

	for(size_t bestGraphNum = 0; bestGraphNum < numBestGraphs; bestGraphNum++) {
		if(bestGraphs[bestGraphNum].getNumEdges() < numEdges) {
			if(nextGamma < 0 || (bestGraphNextGammas[bestGraphNum] > 0 && nextGamma > bestGraphNextGammas[bestGraphNum]))
				nextGamma = bestGraphNextGammas[bestGraphNum];

			double temp = mGamma + (maxScore - bestGraphScores[bestGraphNum]) * 2 / log(mpData->mDataSize) / (numEdges - bestGraphs[bestGraphNum].getNumEdges());
			if(nextGamma < 0 || nextGamma > temp)
				nextGamma = temp;
		}
	}

	if(verbose)
		cout << endl << endl << "< Best Graph >" << endl << bestGraphs.front() << "Number of Edges: " << bestGraphs.front().getNumEdges() << endl << "Next Gamma: " << mNextGamma << endl << endl;

	return(bestGraphs.front());
}
