#include <boost/random.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>


#include "myutil.hpp"


namespace myutil {

	typedef boost::mt19937 generator_t;
	generator_t generator(static_cast<unsigned int>(std::time(0)));
	boost::uniform_int<> uint_dist(0, INT_MAX);
	boost::uniform_real<> ureal_dist(0, 1);
	boost::variate_generator<generator_t&, boost::uniform_int<> > uint_gen(generator, uint_dist);
	boost::variate_generator<generator_t&, boost::uniform_real<> > ureal_gen(generator, ureal_dist);

	size_t rand(size_t range)
	{
		return(uint_gen() % range);
	}

	double rand()
	{
		return(ureal_gen());
	}

	using namespace std;

	bool readDelim(const char* const filename, vector< vector<double> >& data)
	{
		ifstream ifs(filename);
		if(!ifs) {
			cerr << "Cannot read the file. Parse failed." << endl;
			return(false);
		}
		data.clear();

		stringstream iss;
		string line;
		double element;

		while(getline(ifs, line)) {
			if(line.size() == 0)
				break;

			data.resize(data.size() + 1);
			iss << line;
			while(iss >> element)
				data.back().push_back(element);

			if((data.size() > 1) & (data.back().size() != data[data.size() - 2].size())) {
				cerr << "The size of line " << data.size() - 1 << " and line " << data.size() << " are different. Parse failed." << endl;
				return(false);
			}
			iss.clear();
		}

		return(true);
	}

}
