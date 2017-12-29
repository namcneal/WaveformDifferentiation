# include <cstdio>
# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
# include <cassert>
# include <vector>
# include <algorithm>

# include "TCanvas.h"
# include "TROOT.h"
# include "TGraphErrors.h"
# include "TGraph2DErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TH1F.h"
# include "TF1.h"
# include "TLegend.h"
# include "TLatex.h"
# include "TStyle.h"
# include "TApplication.h"
# include "TMultiGraph.h"
# include "TMath.h"
# include "TTree.h"
# include "TFile.h"
# include "TLine.h"
# include "TNtuple.h"
# include "TProfile.h"
# include "TPaveText.h"
# include "TColor.h"

#include "qval_calculation_functions.h"

using namespace std;
const double PERTDC=2.;
const int N_width=15;
const double photon_low=-20.;
const double neutron_low=30.;
const double photon_high=20.;
const double neutron_high=150.;
const double fit_left=-20.;
const double fit_right_long=200.;
const double fit_right=75.;
					  
int main(int argc, char **argv){
	cout << endl << endl << "Loading average waveform data" << endl;
	string dataFileName = argv[2],
		   dataFolderPath = argv[1];
	
	TFile* waveformFile=new TFile((dataFolderPath + dataFileName + "_average_waveforms.root").c_str(), "READ");
	TProfile *avgNeutronWaveform,
			 *avgPhotonWaveform;
	
	waveformFile->GetObject("hist_n_pfx", avgNeutronWaveform);
	waveformFile->GetObject("hist_p_pfx", avgPhotonWaveform);
	
	cout << "Data loaded." << endl;

	double tailStart,
		   tailEnd,
		   neutronQVal,
		   photonQVal;
	
	cout << "Setting analysis parameters" << endl;
	const Int_t numSteps = 10;
	
				
	const double minTailStart = 10,
				 maxTailStart = 80,
				 minTailEnd = 90,
				 maxTailEnd = 200,
				 startStepSize = abs(maxTailStart - minTailStart) / numSteps,
				 endStepSize   = abs(maxTailEnd   - minTailEnd)   / numSteps;
		   
	cout << "Starting Q value calculations" << endl << "\tStep...";
	for(tailStart = minTailStart; tailStart < maxTailStart; tailStart += startStepSize){
		for(tailEnd = minTailEnd; tailEnd < maxTailEnd; tailEnd += endStepSize){

			calculateQVals(tailStart, tailEnd, 
						   avgNeutronWaveform, avgPhotonWaveform,
						   neutronQVal, photonQVal);	
						   
			}
		}
	}
	cout << endl << "Q-value calculations completed." << endl << endl;
	return 0;
}
