#ifndef MYUTIL_HPP_
#define MYUTIL_HPP_

#include <vector>


using namespace std;

typedef unsigned long MyLSizeType;
typedef unsigned short MySizeType;
#define MySizeType_MIN 0
#define MySizeType_MAX USHRT_MAX
#define MyLSizeType_MIN 0
#define MyLSizeType_MAX ULONG_MAX

namespace myutil{
	void myseed(MySizeType);
	MySizeType myrandom(MySizeType);
	double myrandomd();
	double getTSample(MySizeType);
	template<typename T>
	void printvector(vector<T>& vec){
		for(MySizeType i=0 ; i<vec.size() ; i++)
			cout <<vec[i] << " ";
		cout << "\n";
	}
	bool readDelim(const char*, vector< vector<double> >&);
}


#endif /*MYUTIL_HPP_*/
