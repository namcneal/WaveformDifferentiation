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


void calculateQVals(double tailStart, double tailEnd,
					TProfile* avgNeutron, TProfile* avgPhoton,
					double &neutronQVal, double &photonQVal, 
					double &ratioQVal, double &diffQVal);
					  
int main(int argc, char **argv){
	string dataFileName = argv[1],
		   dataFolderPath = "../data/";

	generateAverageWaveforms(dataFileName);
	
	TFile* waveformFile=new TFile((dataFolderPath + dataFileName + "_average_waveforms.root").c_str(), "READ");
	TProfile *avgNeutronWaveform,
			 *avgPhotonWaveform;
			 
	waveformFile->GetObject("hist_n_pfx", avgNeutronWaveform);
	waveformFile->GetObject("hist_p_pfx", avgPhotonWaveform);

	TFile* variedTailFile = new TFile((dataFolderPath + dataFileName + "_varied_tail_qvals.root").c_str(), "RECREATE");
	TTree* qValTree = new TTree("data", "Varied Tail Q-value Results");
	
	double tailStart,
		   neutronQVal,
		   photonQVal,
		   photonNeutronQValRatio,
		   photonNeutronQValDiff;
	   
	qValTree->Branch("tailStart", &tailStart);
	qValTree->Branch("neutronQVal", &neutronQVal);
	qValTree->Branch("photonQVal", &photonQVal);
	qValTree->Branch("ratioQVal", &photonNeutronQValRatio);
	qValTree->Branch("diffQVal", &photonNeutronQValDiff);
	
	const Int_t numSteps = 150,
				minTailStart = 10,
				maxTailStart = 180,
				tailEnd = =170;
	double stepSize = (maxTailStart - minTailStart) / numSteps;
	
	for(tailStart = minTailStart; tailStart <= maxTailStart; tailStart += stepSize){
		
		calculateQVals(tailStart, tailEnd, avgNeutronWaveform, avgPhotonWaveform,
					   neutronQVal, photonQVal, 
					   photonNeutronQValRatio, photonNeutronQValDiff);	
					   
		qValTree->Fill();
	}
	qValTree->Write();
	
	TArrayI params(5);
	params.AddAt(numSteps, 0);	
	params.AddAt(minTailStart, 1);
	params.AddAt(maxTailStart, 2);
	params.AddAt(minTailEnd, 3);
	params.AddAt(maxTailEnd, 4);
	variedTailFile->WriteObjectAny(&params,"TArrayI","calculationParams");
	
	variedTailFile->Close();
	
	return 0;
}


void calculateQVals(double tailStart, double tailEnd, TProfile* avgNeutron, TProfile* avgPhoton,
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
	
	double neutronCurrentQVal = a_tail_n/(a_tail_n+a_main_n),
		   photonCurrentQVal = a_tail_p/(a_tail_p+a_main_p);
		   
	neutronQVal = neutronCurrentQVal;
	photonQVal  = photonCurrentQVal;
							  
}