#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "readdata.hpp"


typedef file_iterator<char> f_iterator_t;
typedef scanner<f_iterator_t> f_scanner_t;
typedef rule<f_scanner_t> f_rule_t;
	
struct TempData
{
    ublas::matrix<double>& data;
    MySizeType colCount;
    MySizeType rowCount;
	vector<double> firstRow;
	vector<MySizeType> firstMissing;
	vector< vector<MySizeType> >& missingIndex;
    TempData(ublas::matrix<double>& data_, vector< vector<MySizeType> >& missingIndex_)
    : data(data_), colCount(0), rowCount(0) , missingIndex(missingIndex_){};
};
     
struct test
{
    test(TempData& temp_data_)
    : temp_data(temp_data_) {}

    void operator()(f_iterator_t first, f_iterator_t last) const
    {
   		   
    }

	TempData& temp_data;
    mutable unsigned i, j;
};

struct update
{
    update(TempData& temp_data_)
    : temp_data(temp_data_) {}

    void operator()(f_iterator_t first, f_iterator_t last) const
    {
    	if(temp_data.rowCount == 0) {
	   		temp_data.data.resize(1000, temp_data.firstRow.size());
	   		temp_data.missingIndex.resize(temp_data.firstRow.size());
	   		for(i = 0; i < temp_data.firstRow.size(); i++)
	   			temp_data.data(0, i) = temp_data.firstRow[i];
	   		for(i = 0; i < temp_data.firstMissing.size(); i++)
	   			temp_data.missingIndex[temp_data.firstMissing[i]].push_back(0);	
    	}
    	else {
	   		if( temp_data.colCount != temp_data.data.size2() ){
	   			cout << "Different number of entries. Parse failed!\n";
	        	exit(-1);
	        } 
	   		if(temp_data.rowCount >= temp_data.data.size1()-1)
	   			temp_data.data.resize(temp_data.data.size1() + 1000, temp_data.data.size2());
    	}
    	temp_data.colCount = 0;
	   	temp_data.rowCount++; 
    }

    TempData& temp_data;
    mutable unsigned i;
};

struct assign_entry
{
    assign_entry(TempData& temp_data_)
    : temp_data(temp_data_) {}

    void operator()(double entry) const
    {
        if(temp_data.rowCount == 0) 
        	temp_data.firstRow.push_back(entry);
        else {
	        temp_data.data(temp_data.rowCount, temp_data.colCount) = entry;
	        if(temp_data.colCount > temp_data.data.size2()){
	        	cout << "Different number of entries. Parse failed!\n";
	        	exit(-1);
	        }
        }
        temp_data.colCount++;
    }

  	TempData& temp_data;
};

struct assign_missing_entry
{
    assign_missing_entry(TempData& temp_data_)
    : temp_data(temp_data_) {}

    void operator()(f_iterator_t first, f_iterator_t last) const
    {
    	if(temp_data.rowCount == 0) {
        	temp_data.firstRow.push_back(0);
        	temp_data.firstMissing.push_back(temp_data.colCount);
    	}
        else {
	        temp_data.data(temp_data.rowCount, temp_data.colCount) = 0;
	        temp_data.missingIndex[temp_data.colCount].push_back(temp_data.rowCount);
	        if(temp_data.colCount > temp_data.data.size2()){
	        	cout << "Different number of entries. Parse failed!\n";
	        	exit(-1);
	        }
        }
        temp_data.colCount++;
    }

  	TempData& temp_data;
};
    
bool readData(const char* filename, ublas::matrix<double>& data, vector< vector<MySizeType> >& missingIndex)
{
	cout << "Calculating covariance matrix...\n";
	
	TempData temp_data(data, missingIndex);
    
    f_iterator_t first(filename);
    if (!first)
    {
        std::cout << "Unable to open file!\n";
        return false;
    }
    f_iterator_t last = first.make_end();


    f_rule_t line = 
    	list_p.direct(real_p[assign_entry(temp_data)]|as_lower_d["na"][assign_missing_entry(temp_data)]
    		, chset<>(" ,\t")) 
    	>> eol_p >> eps_p[update(temp_data)];

    parse_info <f_iterator_t> info = parse(first, last, +line >> *space_p);
    
    if (info.full){
        temp_data.data.resize(temp_data.rowCount, temp_data.data.size2());
        cout << "Parse succeeded!\n";
        cout << temp_data.data.size2() << " variable data of size " << temp_data.data.size1() << "\n\n\n";
        return true;	
    }     
    else{
        std::cout << "Parse failed!\n";
        exit(-1);
    }    
}
