#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "myutil.hpp"
#include "cholesky.hpp"
#include "readdata.hpp" 
#include "graphsampler.hpp"
#include "gjinverse.hpp"


Data::Data(const char* datafile)
{
    	if(!readData(datafile, data, missingIndex)) {
    		cout << "Parse T.resize(data.size2());Failed!\n";
    		exit(-1);
    	}
		nu = 1.0; 
		gamma = 10; 
		lcwishart_cache = vector< vector<double> > (4);    	
    	m = data.size1();
    	n = data.size2();
    	alpha = n;	
    	base_score = -m/2.0*log(2*M_PI) + 1/2.0*log(nu/(m+nu));
    	for(MySizeType i=0 ; i<4 ; i++)
    		lcwishart_cache[i] = vector<double>(n+1, 1);
		arg2[0] = alpha;
		arg2[1] = alpha + m;
		arg2[2] = alpha - 1;
		arg2[3] = alpha + m - 1;
		
		for(MySizeType varNum = 0; varNum < n; varNum++) {
			double sum = 0;
			for(MySizeType rowNum = 0; rowNum < data.size1(); rowNum++)
				sum += data(rowNum, varNum);
			sum /= m - missingIndex[varNum].size();	
			for(vector<MySizeType>::iterator idx = missingIndex[varNum].begin(); idx < missingIndex[varNum].end(); idx++)
				data((*idx), varNum) = sum;
		}			
		
        T0.resize(n, n);
        for(MySizeType rowNum = 0; rowNum < T0.size1(); rowNum++)
        	for(MySizeType colNum = 0; colNum < T0.size2(); colNum++)
        		T0(rowNum, colNum) = (rowNum == colNum) ? 1.0:0.0;
        totalM = ublas::prod(ublas::scalar_matrix<double> (1, (MySizeType)m, 1.0), data);
        totalS = ublas::prod(ublas::trans(data), data);
        T = T0 + totalS - 1 / (nu+m) * ublas::prod(ublas::trans(totalM), totalM);
        test = 0;
};

double Data::getLocalScore(MySizeType node, Structure& strct)
{
	unsigned num_parents  = strct.mParents[node].size();

	if(decompMat.size1() < num_parents + 1)
		decompMat.resize(num_parents + 1, num_parents + 1);

	for(unsigned rowNum = 0; rowNum < num_parents; rowNum++) 
		for(unsigned colNum = 0; colNum <= rowNum; colNum++)
			decompMat(rowNum, colNum) = T(strct.mParents[node][rowNum], strct.mParents[node][colNum]);
	for(unsigned colNum = 0; colNum < num_parents; colNum++)
		decompMat(num_parents, colNum) = T(node, strct.mParents[node][colNum]);
	decompMat(num_parents, num_parents) = T(node, node);	

	if(!cholesky_factorize(decompMat, num_parents + 1)) {
			cout << "Failed Cholesky decomposition!\n";
			exit(-1);
	}

	double val = base_score;
	double det_ratio = pow(decompMat(num_parents, num_parents), 2);
	//assert(det_ratio > 0);

	if(num_parents == 0)
		val += lcwishart(1, 0) - lcwishart(1, 1) - (alpha+m)/2.0*log(det_ratio);
	else{
		val += lcwishart(num_parents+1, 0) + lcwishart(num_parents, 3) - lcwishart(num_parents+1, 1) - lcwishart(num_parents, 2)
				+ (alpha+m)/2.0*(-log(det_ratio)); 
	}
	
	return(val - gamma*num_parents);
}

double Data::lcwishart(MySizeType n, MySizeType index)
{
	if(lcwishart_cache[index][n] < 0)
		return(lcwishart_cache[index][n]);
	
	double val;
	val = n*arg2[index]/2.0*log(2) + n*(n-1)/4*log(M_PI);
	for(MySizeType i=0 ; i<n ; i++)
		val += lgamma((arg2[index]-i)/2.0);
	lcwishart_cache[index][n] = -val;
	return(-val);
}

void Data::imputeData(MySizeType node, Structure& strct)	
{
	vector<MySizeType> mb = strct.getMB(node);
	mb.push_back(node);
	mat mbT(mb.size(), mb.size());
	mat mbMu(1, mb.size());
	double l = m - missingIndex[node].size();
	
	for(MySizeType mbNum = 0; mbNum < mb.size(); mbNum++)
		mbMu(0, mbNum) = totalM(0, mb[mbNum]);
	for(MySizeType indexNum = 0; indexNum < missingIndex[node].size(); indexNum++)
		for(MySizeType mbNum = 0; mbNum < mb.size(); mbNum++)
			mbMu(0, mbNum) -= data(indexNum, mb[mbNum]);	
	totalM(0, node) = mbMu(0, mb.size() - 1);	
	mbMu /= l;		
		
	for(MySizeType mbNum1 = 0; mbNum1 < mb.size(); mbNum1++) {
			for(MySizeType mbNum2 = 0; mbNum2 < mb.size(); mbNum2++) {
				mbT(mbNum1, mbNum2) = totalS(mb[mbNum1], mb[mbNum2]); 
			}
			mbT(mbNum1, mbNum1) += T0(mb[mbNum1], mb[mbNum1]);
	}		
	for(MySizeType indexNum = 0; indexNum < missingIndex[node].size(); indexNum++) {
		for(MySizeType mbNum1 = 0; mbNum1 < mb.size(); mbNum1++) {
			for(MySizeType mbNum2 = 0; mbNum2 < mb.size(); mbNum2++) {
				mbT(mbNum1, mbNum2) -= data(indexNum, mb[mbNum1]) * data(indexNum, mb[mbNum2]); 
			}
		}
		for(MySizeType varNum = 0; varNum < n; varNum++)
			totalS(node, varNum) -= data(indexNum, node) * data(indexNum, varNum);
	}
	mbT -= l*l / (nu+l) * ublas::prod(ublas::trans(mbMu), mbMu);
			
	bool singular;
	mat invT = gjinverse(mbT, singular);
	invT *= (nu+l) / (nu+l+1) * (alpha+l-mb.size()+1);
	assert(!singular);
	mat dataRow(1, mb.size() - 1);
	ublas::matrix_range<mat> invTij(invT, ublas::range(mb.size() - 1, mb.size()), ublas::range(0, mb.size() - 1));
	ublas::matrix_range<mat> invTjj(invT, ublas::range(0, mb.size() - 1), ublas::range(0, mb.size() - 1));
	double invTii = invT(mb.size() - 1, mb.size() - 1);	

	for(MySizeType indexNum = 0; indexNum < missingIndex[node].size(); indexNum++) {
		for(MySizeType mbNum = 0; mbNum < mb.size() - 1; mbNum++)
			dataRow(0, mbNum) = data(indexNum, mb[mbNum]) - mbMu(0, mbNum);
		
		data(indexNum, node) = mbMu(0, mb.size() - 1) 
			- 1 / invTii * (ublas::prod(invTij, ublas::trans(dataRow))) (0,0);

		mat temp = ublas::prod(invTjj - ublas::prod(ublas::trans(invTij), invTij) / invTii, ublas::trans(dataRow));
		double sd = (ublas::prod(dataRow, temp)) (0,0) / invTii;
		test = test > ( sd * invTii / (alpha+l-mb.size() + 1)) ? test : ( sd * invTii / (alpha+l-mb.size() + 1)); 
		sd += (alpha+l-mb.size()+1) / invTii;
		sd /= alpha + l;
		sd = sqrt(sd);
		data(indexNum, node) += sd * myutil::getTSample((MySizeType)(alpha+l));
		
		totalM(0, node) += data(indexNum, node);
		for(MySizeType varNum = 0; varNum < n; varNum++)
			totalS(node, varNum) += data(indexNum, node) * data(indexNum, varNum);
	}

	for(MySizeType varNum = 0; varNum < n; varNum++)
		T(node, varNum) = T0(node, varNum) + totalS(node, varNum) - 1 / (nu+m) *
			totalM(0, varNum) * totalM(0, node);
			
	for(MySizeType mbNum =0; mbNum < mb.size(); mbNum++)
		strct.mLocalScores[mb[mbNum]] = getLocalScore(mb[mbNum], strct);		
}


vector<MySizeType> Structure::getMB(MySizeType node)
{
	vector<MySizeType> mb = mParents[node];
	vector<MySizeType>::iterator new_end;
	
	for(MySizeType childNum = 0; childNum < mChildren[node].size(); childNum++) {
		mb.push_back(mChildren[node][childNum]);
		for(MySizeType spouseNum = 0; spouseNum < mParents[mChildren[node][childNum]].size(); spouseNum++) {
			if(mParents[mChildren[node][childNum]][spouseNum] != node)
				mb.push_back(mParents[mChildren[node][childNum]][spouseNum]);
		}
	}
	
	sort(mb.begin(), mb.end());
	new_end = unique(mb.begin(), mb.end());
	mb.erase(new_end, mb.end());
			
	return(mb);
}	

double Structure::getScore()
{
	double score = 0;
	for(MySizeType i=0 ; i<mNumNodes ; i++)
		score += mLocalScores[i];
	return(score);
}
	
void Structure::randomize(MySizeType num_edges, Data& su){
	MySizeType node1, node2;
	vector<MySizeType>::iterator iter;
	
	for(MySizeType i=0 ; i<mParents.size() ; i++)
		mParents[i].clear();
	for(MySizeType i=0 ; i<mChildren.size() ; i++)
		mChildren[i].clear();	
		
	for(MySizeType i=0 ; i<num_edges ; i++){
		node1 = myutil::myrandom(mNumNodes);
		node2 = myutil::myrandom(mNumNodes-1);
		if(node2 >= node1)
			node2++;
			
		if(node1 < node2){
			if(find(mParents[node1].begin(), mParents[node1].end(), node2) == mParents[node1].end()) {
				mParents[node1].push_back(node2);
				mChildren[node2].push_back(node1);
			}	
		}
		else{	
			if(find(mParents[node2].begin(), mParents[node2].end(), node1) == mParents[node2].end()) {
				mParents[node2].push_back(node1);
				mChildren[node1].push_back(node2);
			}
		}
	}
	
	for(MySizeType i=0; i<mNumNodes; i++)
		mLocalScores[i] = su.getLocalScore(i, *this);
}

bool Structure::isDescendSub(MySizeType node1, MySizeType node2){
	MySizeType i;
	if(isVisited[node1])
		return(0);
	else
		isVisited[node1] = 1;
			
	for(i=0 ; i<mParents[node1].size() ; i++)
		if(mParents[node1][i] == node2)
			return(1);
			
	for(i=0 ; i<mParents[node1].size() ; i++)
		if(isDescendSub(mParents[node1][i], node2))
			return(1);
	
	return(0);
}	
		
bool Structure::isDescendant(MySizeType node1, MySizeType node2){
	MySizeType i;
	
	for(i=0 ; i<mNumNodes ; i++)
		isVisited[i] = 0;
	for(i=0 ; i<mParents[node1].size() ; i++)
		if(mParents[node1][i] != node2)
			if(isDescendSub(mParents[node1][i], node2))
				return(1);
	
	return(0);
}

bool Structure::moveMCMC(Data& su){
	MySizeType node1, node2;
	MySizeType old_edge;
	MySizeType descen;
	bool changed = 0;
	double trans_probs[3];
	double scores[4];
	double rnum;
	double max;
	vector<MySizeType>::iterator iter;
		
	node1 = myutil::myrandom(mNumNodes);
	node2 = myutil::myrandom(mNumNodes-1);
	if(node2 >= node1)
		node2++;
		
	
	if((iter = find(mParents[node1].begin(), mParents[node1].end(), node2)) < mParents[node1].end()){
		old_edge = 1;
		if(isDescendant(node1, node2))
			descen = 1;
		else
			descen = 0;
		
		scores[2] = mLocalScores[node1];
		scores[1] = mLocalScores[node2];
		
		if(descen == 0){
			mParents[node2].push_back(node1);
			scores[3] = su.getLocalScore(node2, *this);
			mParents[node2].pop_back();
		}
		
		mParents[node1].erase(iter);
		scores[0] = su.getLocalScore(node1, *this);	
	}
	else if((iter = find(mParents[node2].begin(), mParents[node2].end(), node1)) < mParents[node2].end()){
		old_edge = 2;
		if(isDescendant(node2, node1))
			descen = 2;
		else
			descen = 0;
		
		scores[3] = mLocalScores[node2];
		scores[0] = mLocalScores[node1];	

		if(descen == 0){
			mParents[node1].push_back(node2);
			scores[2] = su.getLocalScore(node1, *this);
			mParents[node1].pop_back();
		}
			
		mParents[node2].erase(iter);
		scores[1] = su.getLocalScore(node2, *this);
	}
	else{
		old_edge = 0;
		if(isDescendant(node1, node2))
			descen = 1;
		else if(isDescendant(node2, node1))
			descen = 2;
		else
			descen = 0;
				
		scores[0] = mLocalScores[node1];
		scores[1] = mLocalScores[node2];
		
		if(descen != 2){
			mParents[node1].push_back(node2);
			scores[2] = su.getLocalScore(node1, *this);
			mParents[node1].pop_back();
		}
		
		if(descen != 1){
			mParents[node2].push_back(node1);
			scores[3] = su.getLocalScore(node2, *this);
			mParents[node2].pop_back();
		}
	}
	
	if(descen == 0){
		trans_probs[0] = scores[0] + scores[1];
		trans_probs[1] = scores[1] + scores[2];
		trans_probs[2] = scores[0] + scores[3];
	
		max = trans_probs[0];
		if(trans_probs[1] > max)
			max = trans_probs[1];
		if(trans_probs[2] > max)
			max = trans_probs[2];
				
		trans_probs[0] = exp((trans_probs[0] - max)/mTemperature);
		trans_probs[1] = trans_probs[0] + exp((trans_probs[1] - max)/mTemperature);
		trans_probs[2] = trans_probs[1] + exp((trans_probs[2] - max)/mTemperature);
	}
	else if(descen == 1){
		trans_probs[0] = scores[0] + scores[1];
		trans_probs[1] = scores[1] + scores[2];
	
		max = trans_probs[0];
		if(trans_probs[1] > max)
			max = trans_probs[1];
				
		trans_probs[0] = exp((trans_probs[0] - max)/mTemperature);
		trans_probs[1] = trans_probs[0] + exp((trans_probs[1] - max)/mTemperature);
		trans_probs[2] = trans_probs[1];
	}
	else{
		trans_probs[0] = scores[0] + scores[1];
		trans_probs[2] = scores[0] + scores[3];
	
		max = trans_probs[0];
		if(trans_probs[2] > max)
			max = trans_probs[2];
				
		trans_probs[0] = exp((trans_probs[0] - max)/mTemperature);
		trans_probs[1] = trans_probs[0];
		trans_probs[2] = trans_probs[1] + exp((trans_probs[2] - max)/mTemperature);
	}	
	rnum = myutil::myrandomd() * trans_probs[2];
	
	if(rnum < trans_probs[0]){
		if(old_edge == 1) {
			mLocalScores[node1] = scores[0];
			changed = 1;
			mChildren[node2].erase(find(mChildren[node2].begin(), mChildren[node2].end(), node1));
		}
		if(old_edge == 2) {
			mLocalScores[node2] = scores[1];
			changed = 1;
			mChildren[node1].erase(find(mChildren[node1].begin(), mChildren[node1].end(), node2));
		}	
	}
	else if(rnum < trans_probs[1]){
		if(old_edge == 0) {
			mLocalScores[node1] = scores[2];
			mChildren[node2].push_back(node1);
			changed = 1;
		}
		if(old_edge == 2) {
			mLocalScores[node1] = scores[2];
			mLocalScores[node2] = scores[1];
			changed = 1;
			mChildren[node1].erase(find(mChildren[node1].begin(), mChildren[node1].end(), node2));
			mChildren[node2].push_back(node1);
		}	
		mParents[node1].push_back(node2);
	}
	else{
		if(old_edge == 0) {
			mLocalScores[node2] = scores[3];
			mChildren[node1].push_back(node2);
			changed = 1;
		}
		if(old_edge == 1) {
			mLocalScores[node1] = scores[0];
			mLocalScores[node2] = scores[3];
			changed = 1;
			mChildren[node2].erase(find(mChildren[node2].begin(), mChildren[node2].end(), node1));
			mChildren[node1].push_back(node2);
		}
		mParents[node2].push_back(node1);
	}
	
	return(changed);			
}