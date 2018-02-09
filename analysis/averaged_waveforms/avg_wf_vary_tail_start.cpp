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

#include "../qval_calculation_functions.h"

int main(int argc, char **argv){
	cout << endl << endl << "Loading average waveform data" << endl;
	string dataFileName = argv[2],
		   dataFolderPath = argv[1];
	
	TFile* waveformFile=new TFile((dataFolderPath + "/averaged_waveforms/" + dataFileName + "_average_waveforms.root").c_str(), "READ");
	if(waveformFile->IsZombie()){
		cout << endl << "Could not read the averaged waveform file" << endl;
		return -1;
	}

	TProfile *avgNeutronWaveform,
			 *avgPhotonWaveform;
	
	waveformFile->GetObject("hist_n_pfx", avgNeutronWaveform);
	waveformFile->GetObject("hist_p_pfx", avgPhotonWaveform);
	cout << "Data loaded." << endl;


	cout << "Setting analysis parameters" << endl;
	const int NUM_STEPS = 10,
			  MIN_TAIL_START = 10,
			  MAX_TAIL_START = 80;

	const double START_STEP_SIZE = abs(MIN_TAIL_START - MAX_TAIL_START ) / NUM_STEPS;
	
 	cout << "Creating the file to hold the Q-value data" << endl;
	double tailStart = 0,
		   tailEnd   = 170,
		   neutronQVal = 0,
		   photonQVal  = 0;

	stringstream ss;
	ss << dataFolderPath << "/averaged_waveforms/" << "qvalues_varying_start_from_"
	   << MIN_TAIL_START << "_to_" << MAX_TAIL_START << "_with_tail_end_at_" << tailEnd << ".root";
	
	string fileName = ss.str();

	TFile  * qValFile = new TFile(fileName.c_str(), "RECREATE");
	TTree *tree = new TTree("Averaged Waveform Q-Value Data Varying Tail", "Q-Value Data");

	tree->Branch("Start of Tail",  &tailStart,   "tailStart/D");
	tree->Branch("Neuton Q-Value", &neutronQVal, "neutronQVal/D");
	tree->Branch("Photon Q-Value", &photonQVal,   "photonQVal/D");
		   
	cout << "Starting Q value calculations" << endl << "\tStep...";
	for(tailStart = MIN_TAIL_START; tailStart < MAX_TAIL_START; tailStart += START_STEP_SIZE){
			
			calculateQVal(tailStart, tailEnd, 
						   avgNeutronWaveform, avgPhotonWaveform,
						   neutronQVal, photonQVal);	

			tree->Fill();
						   
	}
	cout << endl << "Q-value calculations completed." << endl << "Writing to file." << endl;

	tree->Write();

	waveformFile->Close();
	qValFile->Close();
	return 0;
}
