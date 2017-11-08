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

using namespace std;
const double PERTDC=2.;
const int N_width=15;

int main(int argc, char **argv){
	gROOT->Reset();
	gStyle->SetOptStat(110010);
	gErrorIgnoreLevel = kWarning;
	string name (argv[1]);
	
	string plot_name=name+"_spectrum.pdf";
	TCanvas* c0=new TCanvas("0", "", 800, 600);
	
	TFile* raw_file=new TFile((name+"_raw_data.root").c_str(), "READ");
	TTree* raw_tree=(TTree*)raw_file->Get(Form("%s_raw", name.c_str()));
	
	double ped_n, amp_n;
	int is_good_n;
	int adc_n[N_width*32];
	int n_bad_data=0;
	
	raw_tree->SetBranchAddress("adc_n", adc_n);
	raw_tree->SetBranchAddress("n_quality", &is_good_n);
	
	int n_tot=raw_tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<name<<"_raw_data.root"<<endl;
	
	TH1D* h_spec=new TH1D("h_spec", "Spectrum; ADC; count", 200, 0., 4000.);
	
	for (int i=0; i<n_tot; i++){
		ped_n=-999.;
		amp_n=-999.;
		
		raw_tree->GetEntry(i);
		
		if (is_good_n==0) n_bad_data++;
		else{
			ped_n=0.;
			for (int i=0; i<10; i++) ped_n+=(double)adc_n[i];
			ped_n/=10.;
			
			for (int i=0; i<400; i++){
				if ((double)adc_n[i]-ped_n>(amp_n)) {
					amp_n=(double)adc_n[i]-ped_n;
				}
			}
			h_spec->Fill(amp_n);
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_bad_data<<" events have bad data."<<endl;
	
	gStyle->SetOptStat(110010);
	h_spec->SetLineColor(kRed);
	
	h_spec->Draw();
	c0->Print((plot_name+"(").c_str(),"pdf");
	
	c0->SetLogy();
	c0->Print((plot_name).c_str(),"pdf");
	/* Na-22
	h_spec->DrawClone();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TF1 *f1=new TF1("f1", "[0]+[1]*TMath::Gaus(x, [2], [3], 1)", 300, 550);
	f1->SetParameters(0., 1000., 350., 20.);
	f1->SetParLimits(3, 1., 200.);
	f1->SetParNames("A_{0}", "A", "#mu", "#sigma");
	TFitResultPtr frp1=h_spec->Fit(f1, "QRSL");
	f1->SetLineColor(kBlue);
	f1->DrawClone("SAME");
	c0->Print((plot_name).c_str(),"pdf");
	
	h_spec->DrawClone();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TF1 *f2=new TF1("f2", "[0]+[1]*TMath::Gaus(x, [2], [3], 1)", 1200, 2000);
	f2->SetParameters(0., 1000., 1300., 200.);
	f2->SetParNames("A_{0}", "A", "#mu", "#sigma");
	TFitResultPtr frp2=h_spec->Fit(f2, "QRSL");
	f2->SetLineColor(kBlue);
	f2->DrawClone("SAME");
	c0->Print((plot_name+")").c_str(),"pdf");
	*/
	/* Cs-137
	h_spec->DrawClone();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TF1 *f2=new TF1("f2", "[0]+[1]*TMath::Gaus(x, [2], [3], 1)", 400, 800);
	f2->SetParameters(0., 1000., 400., 20.);
	f2->SetParLimits(3, 1., 200.);
	f2->SetParNames("A_{0}", "A", "#mu", "#sigma");
	TFitResultPtr frp2=h_spec->Fit(f2, "QRSL");
	f2->SetLineColor(kBlue);
	f2->DrawClone("SAME");
	c0->Print((plot_name+")").c_str(),"pdf");
	*/
	// /* Co-60
	h_spec->DrawClone();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	TF1 *f2=new TF1("f2", "[0]+[1]*TMath::Gaus(x, [2], [3], 1)", 1100, 2000);
	f2->SetParameters(0., 1000., 1100., 200.);
	f2->SetParNames("A_{0}", "A", "#mu", "#sigma");
	TFitResultPtr frp2=h_spec->Fit(f2, "QRSL");
	f2->SetLineColor(kBlue);
	f2->DrawClone("SAME");
	c0->Print((plot_name+")").c_str(),"pdf");
	// */
	return 0;
}


