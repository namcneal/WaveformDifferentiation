#include <iostream>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "../qval_calculation_functions.h"

using namespace std;
void qvals_vary_start_with_end_at_170(string dataFolderPath, string dataFileName){
	
	string avgWaveformROOTFileName = dataFolderPath + "/" + dataFileName + "_average_waveforms.root";

	TFile * waveformFile = new TFile(avgWaveformROOTFileName.c_str(), "READ");
		
	if(!waveformFile){
		cout << "Could not find/read the waveform file" << endl;
	}

	cout << "Read in the waveform file" << endl;

	TProfile * avgNeutronWaveform = NULL,
			 * avgPhotonWaveform  = NULL;
	waveformFile->GetObject("hist_n_pfx", avgNeutronWaveform);
	waveformFile->GetObject("hist_p_pfx", avgPhotonWaveform);

	int numSteps = 100,
		minTailStart = 25,
		maxTailStart  = 120,
		tailEnd      = 170;

	double stepSize = (maxTailStart - minTailStart) / (double)numSteps;

	double tailStart = 0,
		   neutronQVal = 0,
		   photonQVal  = 0;

	TGraph  *neutronGraph = new TGraph(),
			*photonGraph  = new TGraph();

	int i = 0;
	for(tailStart = minTailStart; tailStart <= maxTailStart; tailStart += stepSize){
		
		calculateQVal(tailStart, tailEnd, 
					  avgNeutronWaveform, avgPhotonWaveform,
					  neutronQVal, photonQVal);
		cout << photonQVal << endl;
		neutronGraph->SetPoint(i, tailStart, neutronQVal);
		photonGraph ->SetPoint(i, tailStart, photonQVal );
		
		i++;
	}
	
	TCanvas *c0 = new TCanvas("c0", "c0",800,500);
	TMultiGraph *mg = new TMultiGraph();
	
	neutronGraph->SetLineColor(kRed);
	mg->Add(neutronGraph, "l");
	
	photonGraph->SetLineColor(kBlue);
	mg->Add(photonGraph,"l");
	
	mg->Draw("a");
	
	 c0->Modified();
	 c0->Update();
//	
}