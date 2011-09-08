#include "graphsampler.hpp"
#include "readdata.hpp"
#include "myutil.hpp"

#include <fstream>

int main(int argc, char *argv[]){
	string datafile;
	MySizeType num_chains, j, k, swap, temp1, swap_period, impute_period;
	MyLSizeType i, burnin, num_iter;
    double energy, temp, max_temp;
    fstream ouf;
    
    if(argc < 2){
    	cout << "Datafile: ";
    	cin >> datafile;
    }
    else 
    	datafile = string(argv[1]);
    Data su(datafile.c_str());
    
    
    cout << "Gamma: ";
    cin >> su.gamma;
    cout << "Number of chains: ";
    cin >> num_chains;
    cout << "Maximum temperature: ";
    cin >> max_temp;
    cout << "Swap period: ";
    cin >> swap_period;
    cout << "Imputation period: ";
    cin >> impute_period;
    cout << "Burnin (k): ";
	cin >> burnin;
	burnin *= 1000;
	cout << "Number of iterations (k): ";
	cin >> num_iter;
	num_iter *= 1000;

    vector<double> swap_rate(num_chains-1);
    vector<double> mcmc_rate(num_chains);
	vector<unsigned short> idx(num_chains);

	myutil::myseed(1574);
	

    for(j=0 ; j<num_chains ; j++)
    	idx[j] = j;
    	
    Structure strct[num_chains];
    for(j=0 ; j<num_chains ; j++)
    	strct[j] = Structure(su);
    
    vector< vector<double> > edges;
    for(MySizeType j=0 ; j<su.getNumNodes() ; j++)
    	edges.push_back(vector<double>(su.getNumNodes()));

    strct[0].mTemperature = 1.;
    for(j=1 ; j<num_chains ; j++)
    	strct[j].mTemperature = strct[j-1].mTemperature * pow(max_temp, 1./(num_chains-1));


	cout << "\nRunning chains with temperatures: ";
	for(j=0 ; j<num_chains ; j++)
    	cout << strct[j].mTemperature << " ";
    cout << "\n\n";	
	
	ouf.open("output.txt", ios::out);
	if(ouf == NULL){
		cout << "Failed to write\n";
		exit(-1);
	}

    
    
    cout << "LScore: ";
    printf("%.10e\n", strct[idx[0]].getScore()); 
    for(j=0 ; j<su.getNumNodes() ; j++){
	    	cout << j << ": ";
	    	myutil::printvector(strct[0].mParents[j]);
    }
 	cout << "\n\n";
 	   
    for(i=1 ; i<=burnin ; i++){
		if(i%10000 == 0){
			cout << "burnin: " << i/1000 << "K\n";
			cout << "LScore: "; 
			printf("%.10e\n", strct[idx[0]].getScore());
			for(j=0 ; j<swap_rate.size() ; j++)
				swap_rate[j] /= (10000/swap_period);
			for(j=0 ; j<mcmc_rate.size() ; j++)
				mcmc_rate[j] /= 10000;
			cout << "swap rates: ";
			myutil::printvector(swap_rate);
			cout << "mcmc rates: ";
			for(j=0 ; j<num_chains ; j++)
				cout << mcmc_rate[idx[j]] << " ";
			cout << "\n";	
			for(j=0 ; j<swap_rate.size() ; j++)
				swap_rate[j] = 0;
			for(j=0 ; j<mcmc_rate.size() ; j++)
				mcmc_rate[j] = 0;
			cout << "max ratio: " << su.test;	
			cout << "\n\n";
			ouf.close();
			ouf.open("output.txt", ios::app | ios::out);
		}
		if(i%impute_period == 0)
			su.imputeData(myutil::myrandom(su.getNumNodes()), strct[idx[0]]);
		if(i%swap_period == 0){
			for(j=0 ; j<num_chains-1 ; j++){
				swap = num_chains-2-j;
				energy = (1/strct[idx[swap]].mTemperature - 1/strct[idx[swap+1]].mTemperature)*(-strct[idx[swap]].getScore() + strct[idx[swap+1]].getScore());
				if(energy > 0)
					energy = 0;
				energy = exp(energy);
				if(myutil::myrandomd() < energy){
					swap_rate[swap]++;
					temp = strct[idx[swap]].mTemperature;
					strct[idx[swap]].mTemperature = strct[idx[swap+1]].mTemperature;
					strct[idx[swap+1]].mTemperature = temp; 
					temp1 = idx[swap];
					idx[swap] = idx[swap+1];
					idx[swap+1] = temp1;
				}	
			}
		}	
	    
	    ouf << 	strct[idx[0]].getScore();	
	    for(j=1 ; j<num_chains ; j++)
	    	ouf <<"\t" << strct[idx[j]].getScore();
	    for(j=0 ; j<num_chains ; j++){
    		if(strct[idx[j]].moveMCMC(su))
    			mcmc_rate[idx[j]]++;
	    }
	    ouf << "\n";
    }

    for(i=1 ; i<=num_iter ; i++){
		if(i%10000 == 0){
			cout << "iter: " << i/1000 << "K\n";
			cout << "LScore: ";
			printf("%.10e\n", strct[idx[0]].getScore());
			for(j=0 ; j<swap_rate.size() ; j++)
				swap_rate[j] /= (10000/swap_period);
			for(j=0 ; j<mcmc_rate.size() ; j++)
				mcmc_rate[j] /= 10000;
			cout << "swap rates: ";
			myutil::printvector(swap_rate);
			cout << "mcmc rates: ";
			for(j=0 ; j<num_chains ; j++)
				cout << mcmc_rate[idx[j]] << " ";
			cout << "\n";	
			for(j=0 ; j<swap_rate.size() ; j++)
				swap_rate[j] = 0;
			for(j=0 ; j<mcmc_rate.size() ; j++)
				mcmc_rate[j] = 0;
			cout << "max ratio: " << su.test;	
			cout << "\n\n";
			
			ouf.close();
			ouf.open("output.txt", ios::app | ios::out);
			
			std::ofstream ouf1("edges.txt");
			if(ouf1 == NULL){
					cout << "Failed to write\n";
					exit(-1);
			}
			for(j=0 ; j<edges.size() ; j++){
				ouf1 << edges[j][0] / i;
				for(k=1 ; k<edges.size() ; k++)
					ouf1 << "\t" << edges[j][k] / i;
				ouf1 << "\n";
			}
			ouf1.close();
		}
		if(i%impute_period == 0)
			su.imputeData(myutil::myrandom(su.getNumNodes()), strct[idx[0]]);
		if(i%swap_period == 0){
			for(j=0 ; j<num_chains-1 ; j++){
				swap = num_chains-2-j;
				energy = (1/strct[idx[swap]].mTemperature - 1/strct[idx[swap+1]].mTemperature)*(-strct[idx[swap]].getScore() + strct[idx[swap+1]].getScore());
				if(energy > 0)
					energy = 0;
				energy = exp(energy);
				if(myutil::myrandomd() < energy){
					swap_rate[swap]++;
					temp = strct[idx[swap]].mTemperature;
					strct[idx[swap]].mTemperature = strct[idx[swap+1]].mTemperature;
					strct[idx[swap+1]].mTemperature = temp; 
					temp1 = idx[swap];
					idx[swap] = idx[swap+1];
					idx[swap+1] = temp1;
				}	
			}
		}		
			
	    for(j=0 ; j<num_chains ; j++){
    		if(strct[idx[j]].moveMCMC(su))
    			mcmc_rate[idx[j]]++;
	    }
	    
	    for(j=0 ; j<su.getNumNodes() ; j++) {
	    	for(k=0 ; k<strct[idx[0]].mParents[j].size() ; k++) {
	    		edges[j][strct[idx[0]].mParents[j][k]]++;
	    		edges[strct[idx[0]].mParents[j][k]][j]++;
		}
	    }	
	    		
		ouf << 	(strct[idx[0]].getScore());	
	    for(j=1 ; j<num_chains ; j++)
	    	ouf <<"\t" << strct[idx[j]].getScore();
	    ouf << "\n";	
    }
    
    
    cout << "\nLScore: ";
    printf("%.10e\n", strct[idx[0]].getScore());
    for( j=0 ; j<su.getNumNodes() ; j++){
	    	cout << j << ": ";
	    	myutil::printvector(strct[idx[0]].mParents[j]);
    }    

	ouf.close();
	
}

