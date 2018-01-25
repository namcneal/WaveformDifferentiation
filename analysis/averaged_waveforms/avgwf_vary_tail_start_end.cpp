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

void calculateQVal(double tailStart, double tailEnd, 
				   TProfile* avgNeutron, TProfile* avgPhoton,
				   double &neutronQVal, double &photonQVal);
					  
void avgwf_vary_tail_start_end(string dataFolderPath, string dataFileName){	
	string full_file_path = dataFolderPath + '/' + dataFileName + "_average_waveforms.root";
	cout << endl << endl <<"Loading averaged waveform data from " << full_file_path << endl;
	
	TFile* waveformFile = TFile::Open(full_file_path.c_str());
	if(!waveformFile){
		cout << "Could not open " << full_file_path << endl;
	}
	
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
				 
	TGraph2D *neutronQValGraph = new TGraph2D(numSteps),
			 *photonQValGraph  = new TGraph2D(numSteps);
	neutronQValGraph->SetName("Neutron Q-Values");
	photonQValGraph->SetName("Photon Q-Values");
	
	int i = 0;
	cout << "Starting Q value calculations" << endl << "\tStep...";
	for(tailStart = minTailStart; tailStart < maxTailStart; tailStart += startStepSize){
		for(tailEnd = minTailEnd; tailEnd < maxTailEnd; tailEnd += endStepSize){

			calculateQVal(tailStart, tailEnd, 
						  avgNeutronWaveform, avgPhotonWaveform,
						  neutronQVal, photonQVal);

			neutronQValGraph->SetPoint(i, tailStart, tailEnd, neutronQVal);
			 photonQValGraph->SetPoint(i, tailStart, tailEnd, photonQVal );
									
		}
		
		i++;

	}
	cout << endl << "Q-value calculations completed." << endl << endl;
	
	TCanvas *diffCanvas = new TCanvas("qValDiff", "Q-value Relative Difference", 700, 600);
	neutronQValGraph->Draw("COLZ");
	

	diffCanvas->SaveAs((dataFolderPath + dataFileName + "_varied_start_tail_difference.png").c_str(), "png");
}

void calculateQVal(double tailStart, double tailEnd, 
				   TProfile* avgNeutron, TProfile* avgPhoton,
				   double &neutronQVal, double &photonQVal){
	
	int tailStartBin = avgNeutron->FindBin(tailStart),
		tailEndBin   = avgNeutron->FindBin(tailEnd);
		
	double a_main_n=0.,
		   a_tail_n=0.,
		   a_main_p=0.,
		   a_tail_p=0.;
	
	for (int i=1; i<tailStartBin; i++){
		a_main_n+=avgNeutron->GetBinContent(i);
		a_main_p+=avgPhoton->GetBinContent(i);
	}
	for (int i = tailEndBin; i<=avgNeutron->FindBin(199.75); i++){
		a_tail_n+=avgNeutron->GetBinContent(i);
		a_tail_p+=avgPhoton->GetBinContent(i);
	}
	
	double neutronCurrentQVal = a_tail_n / (a_tail_n+a_main_n),
		   photonCurrentQVal  = a_tail_p / (a_tail_p+a_main_p);
		   
	neutronQVal = neutronCurrentQVal;
	photonQVal  = photonCurrentQVal;
}
