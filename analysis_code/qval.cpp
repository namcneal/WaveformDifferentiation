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
const double tail_start=20.;
const double fit_left=-20.;
const double fit_right_long=200.;
const double fit_right=75.;

int main(int argc, char **argv){	
	gStyle->SetOptStat(111111);
	gErrorIgnoreLevel = kWarning;

	string dataFolderPrefix = "../data/",
		   dataFileName = argv[1];
	
	TFile* par_file=new TFile((dataFolderPrefix + dataFileName + "_coarse_time.root").c_str(), "UPDATE");
	TTree* par_tree=(TTree*)par_file->Get((dataFileName + "_coarse_time").c_str());
	if (!par_tree){
		cout << "Could not find the " <<  dataFileName << "_coarse_time.root file" << endl;
	}
	
	TFile* raw_file=new TFile((dataFolderPrefix + dataFileName + "_raw_data.root").c_str(), "READ");
	TTree* raw_tree=(TTree*)raw_file->Get(Form("%s_raw", dataFileName.c_str()));
	if (!raw_tree){
			cout << "Could not find the " <<  dataFileName << "raw_time.root file" << endl;
		}

	
	double ped_n, ped_f, amp_n, amp_f, t_max_n, t_max_f, t_50_n, t_50_f;
	int is_good_n = NULL, 
		is_good_f = NULL,
		adc_n[N_width*32] = {NULL},
		adc_f[N_width*32] = {NULL};
	
	cout << "Setting parse tree branch addresses" << endl;
	par_tree->SetBranchAddress("ped_n", &ped_n);
	par_tree->SetBranchAddress("ped_f", &ped_f);
	par_tree->SetBranchAddress("amp_n", &amp_n);
	par_tree->SetBranchAddress("amp_f", &amp_f);
	par_tree->SetBranchAddress("t_max_n", &t_max_n);
	par_tree->SetBranchAddress("t_max_f", &t_max_f);
	par_tree->SetBranchAddress("t_50_n", &t_50_n);
	par_tree->SetBranchAddress("t_50_f", &t_50_f);

	cout << "Setting raw tree branch addresses" << endl;
	raw_tree->SetBranchAddress("adc_n", adc_n);
	raw_tree->SetBranchAddress("adc_f", adc_f);
	raw_tree->SetBranchAddress("n_quality", &is_good_n);
	raw_tree->SetBranchAddress("f_quality", &is_good_f);
	

	cout << "Getting entries.";
	int n_tot=raw_tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<dataFileName<<"_raw_data.root"<<endl;

	TH1D* qval=new TH1D("qval","Distribution of tail Q/total Q; tail Q/total Q; counts", 100, 0., 1);
	TH1D* qval_n=new TH1D("qval_n","Distribution of tail Q/total Q; tail Q/total Q; counts", 100, 0., 1);
	TH1D* qval_p=new TH1D("qval_p","Distribution of tail Q/total Q; tail Q/total Q; counts", 100, 0., 1);

	int n_neutron=0;
	int n_photon=0;

	cout << "Starting loop";
	for (int i=0; i<n_tot; i++){
		par_tree->GetEntry(i);
		raw_tree->GetEntry(i);

		if (is_good_f==1 && is_good_n==1 && t_max_f!=0. && t_max_n!=0.){
			double t_50_diff=t_50_f-t_50_n;
			double start_time=t_50_f-20.;
			int start_TDC=(int)ceil(start_time);
			if (start_TDC%2) start_TDC+=1;
			start_TDC/=2;

			int tail_TDC=start_TDC+20;
			int end_TDC=start_TDC+110;

			double a_all=0.;
			double a_tail=0.;

			for (int j=start_TDC; j<tail_TDC; j++){
				a_all+=adc_f[j]-ped_f;
			}

			// The integration of the waveforms.
			for (int j=tail_TDC; j<end_TDC; j++){
				a_all+=(double)adc_f[j]-ped_f;
				a_tail+=(double)adc_f[j]-ped_f;
			}
			double ratio=a_tail/a_all;
			qval->Fill(ratio);

			if (t_50_diff>neutron_low && t_50_diff<neutron_high){
				qval_n->Fill(ratio);
				n_neutron++;
			}
			if (t_50_diff>photon_low && t_50_diff<photon_high){
				qval_p->Fill(ratio);
				n_photon++;
			}
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_neutron<<" neutron events are present."<<endl;
	cout<<"    "<<n_photon<<" photon events are present."<<endl;
	
	const Int_t NRGBs = 5;
	const Int_t NCont = 256;
	Double_t stops[NRGBs] = { 0.00, 0.30, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs] = { 0.00, 0.00, 0.57, 0.90, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.65, 0.95, 0.20, 0.00 };
	Double_t blue[NRGBs] = { 0.51, 0.55, 0.15, 0.00, 0.10 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
	
	string plot_name = dataFolderPrefix + dataFileName + "_qval.pdf";
	TCanvas* c0=new TCanvas("0", "", 1600, 1200);
	c0->Print((plot_name+"[").c_str(),"pdf");
	
	qval->SetMarkerColor(kBlack);
	qval_n->SetMarkerColor(kRed);
	qval_p->SetMarkerColor(kBlue);
	qval->SetLineColor(kBlack);
	qval_n->SetLineColor(kRed);
	qval_p->SetLineColor(kBlue);

	qval_n->Draw();
	TLegend* leg1=new TLegend(.75, .60, .9, .55, "");
	leg1->SetFillColor(0);
	leg1->AddEntry(qval_n, "neutron");
	leg1->Draw("SAMES");
	c0->Print(plot_name.c_str(),"pdf");

	qval_p->Draw();
	TLegend* leg2=new TLegend(.75, .60, .9, .55, "");
	leg2->SetFillColor(0);
	leg2->AddEntry(qval_p, "photon");
	leg2->Draw("SAMES");
	c0->Print(plot_name.c_str(),"pdf");

	TCanvas* c1=new TCanvas("1", "", 1600, 1200);
	gStyle->SetOptStat(0);
	qval->Draw();
	qval_n->Draw("SAME");
	qval_p->Draw("SAME");
	TLegend* leg3=new TLegend(.75, .650, 1, .95, "");
	leg3->SetFillColor(0);
	leg3->AddEntry(qval_p, "photon");
	leg3->AddEntry(qval_n, "neutron");
	leg3->AddEntry(qval, "sum");
	leg3->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");

	c1->SetLogy();
	c1->Print((plot_name+")").c_str(),"pdf");

	par_file->Close();
	raw_file->Close();
	return 0;
}