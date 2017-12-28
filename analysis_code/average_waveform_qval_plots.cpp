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

int main(int argc, char **argv){
	cout << "Reading in the data and creating the trees" << endl;
	string dataFileName = argv[1],
		   dataFolderPath = "../data/";
		   
	string plot_name = dataFolderPath + dataFileName + "PLEASE.pdf";
	
	TFile* readFile=new TFile((dataFolderPath + dataFileName+"_varied_tail_qvals.root").c_str(), "READ");
	TTree* readTree=(TTree*)readFile->Get("data");

	TProfile  *avgNeutron =  0,
			  *avgPhoton  = 0;

	TArrayI calculationParams(5);

	readFile->GetObject("hist_n_pfx", avgNeutron);
	readFile->GetObject("hist_p_pfx", avgPhoton);

	readFile->GetListOfKeys()->Print();
	readFile->GetObject("calculationParams", calculationParams);


	// double tailStarts[numSteps], 
	// 	   neutronQVals[numSteps],
	// 	   photonQVals[numSteps],
	// 	   photonNeutronQValRatios[numSteps],
	// 	   photonNeutronQValDiffs[numSteps];
	
	
	// for(int i =0; i <= numSteps; i++){
	// 	readTree->SetBranchAddress("tailStart", tailStarts + i);
	// 	readTree->SetBranchAddress("neutronQVal", neutronQVals + i);
	// 	readTree->SetBranchAddress("photonQVal", photonQVals + i);
	// 	readTree->SetBranchAddress("ratioQVal", photonNeutronQValRatios + i);
	// 	readTree->SetBranchAddress("diffQVal", photonNeutronQValDiffs + i);
		
	// 	readTree->GetEntry(i);
	// }
	
	// TCanvas* c0=new TCanvas("0", "", 600, 500);
	// c0->Divide(1,2);
	
	// c0->cd(1);
	// TGraph* neutronPlot = new TGraph(numSteps, tailStarts, neutronQVals);
	// neutronPlot->SetTitle("Photon (Blue) & Neutron (Red) Individual Q-values and Difference");
	// neutronPlot->GetXaxis()->SetLimits(0, 200);
	// neutronPlot->SetMinimum(0);
	// neutronPlot->SetMaximum(.6);
	// neutronPlot->SetMarkerStyle(24);
	// neutronPlot->SetMarkerSize(.5);
	// neutronPlot->SetMarkerColor(2);
	// neutronPlot->SetLineColor(2);
	// neutronPlot->Draw("Al");
	
	// TGraph* photonPlot = new TGraph(numSteps, tailStarts, photonQVals);
	// photonPlot->SetMarkerStyle(24);
	// photonPlot->SetMarkerColor(4);
	// photonPlot->SetMarkerSize(.5);
	// photonPlot->SetLineColor(4);
	// photonPlot->Draw("l");
	
	// TGraph* ratioPlot = new TGraph(numSteps, tailStarts, photonNeutronQValDiffs);
	// ratioPlot->SetMarkerStyle(24);
	// ratioPlot->SetMarkerColor(3);
	// ratioPlot->SetMarkerSize(.5);
	// ratioPlot->SetLineColor(3);
	// ratioPlot->Draw("l");
	
	// TLine* line_top =new TLine(20,0,20,1.1);
	// line_top->SetLineColor(kBlack);
	// line_top->Draw("SAME");
		
	// c0->cd(2);
	// avgNeutron->SetLineColor(kRed);
	// avgPhoton->SetLineColor(kBlue);
	// avgNeutron->SetLineWidth(3);
	// avgPhoton->SetLineWidth(3);
	// avgNeutron->SetStats(false);
	// avgPhoton->SetStats(false);
	// avgNeutron->GetXaxis()->SetLimits(0, 200);
	// avgPhoton->GetXaxis()->SetLimits(0, 200);
	// avgNeutron->Draw("HIST C");
	// avgPhoton->Draw("A HIST SAME C");
	
	// TLine* line_bottom =new TLine(20,0,20,1.1);
	// line_bottom->SetLineColor(kBlack);
	// line_bottom->Draw("SAME");
	
	
	// c0->Modified();
	// c0->Update();
	// c0->SaveAs((dataFolderPath + dataFileName + "_varied_qvals_DIFFERENCE.png").c_str(), "png");
	
	readFile->Close();
	return 0;
}