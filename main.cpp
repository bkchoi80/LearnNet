#include <fstream>
#include <sstream>
#include <list>

#include "myutil.hpp"
#include "learnnet.hpp"


using namespace std;

int main(int argc, char *argv[]){
	string datafile;
	string outputfile;
	vector< vector<double > > data;
	if(argc < 2){
    	cout << "Datafile: ";
    	cin >> datafile;
    }
    else
    	datafile = string(argv[1]);
	while(!myutil::readDelim(datafile.c_str(), data)) {
		cout << "Cannot read \"" << datafile << "\". Please reenter: ";
		cin >> datafile;
	}

	size_t numIter, tabuPeriod, restartPeriod;
	double dataSize, gamma;

	cout << "Number of Data (k): ";
    cin >> dataSize;
    dataSize *= 1000;
    cout << "Initial Gamma: ";
    cin >> gamma;
    cout << "Tabu period: ";
    cin >> tabuPeriod;
    cout << "Restart period (k): ";
    cin >> restartPeriod;
    restartPeriod *= 1000;
    cout << "Number of iterations (k): ";
	cin >> numIter;
	numIter *= 1000;
	cout << "Output file: ";
	cin >> outputfile;

	std::ofstream ouf;
	Data inputData(data, dataSize);
	OrderSearch orderSearch(&inputData, gamma);
	orderSearch.randomize();

	do {
		Dag bestGraph = orderSearch.doTabuSearch(numIter, tabuPeriod, restartPeriod, true);

		ouf.open(outputfile.c_str(), ios_base::app);
		ouf << orderSearch.getGamma() << endl;
		ouf << bestGraph;
		ouf.close();

		if(bestGraph.getNumEdges() == 0)
			break;

		orderSearch.moveToNextGamma();
	} while(1);
}

