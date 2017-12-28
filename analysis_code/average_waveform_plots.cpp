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

	TFile* readFile=new TFile((dataFolderPath + dataFileName+"_average_waveforms.root").c_str(), "READ");
	TTree* readTree=(TTree*)readFile->Get("data");

	TProfile  *avgNeutron =  0,
			  *avgPhoton  = 0;

	TArrayI calculationParams(5);

	readFile->GetObject("hist_n_pfx", avgNeutron);
	readFile->GetObject("hist_p_pfx", avgPhoton);

	 TCanvas* c0=new TCanvas("0", "", 1500, 1000);
	 gStyle->SetOptTitle(0);
	 
	 avgNeutron->SetLineColor(kRed);
	 avgPhoton->SetLineColor(kBlue);
	 avgNeutron->SetLineWidth(3);
	 avgPhoton->SetLineWidth(3);
	 avgNeutron->SetStats(false);
	 avgPhoton->SetStats(false);
	 avgNeutron->GetXaxis()->SetLimits(0, 200);
	 avgPhoton->GetXaxis()->SetLimits(0, 200);
	 avgNeutron->Draw("HIST C");
	 avgPhoton->Draw("A HIST SAME C");

	 c0->Modified();
	 c0->Update();
	 c0->SaveAs((dataFolderPath + dataFileName + "avg_waveforms.png").c_str(), "png");

	readFile->Close();
	return 0;
}
