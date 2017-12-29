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

const string DATA_FOLDER_PATH = "../data/";
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

/*
Requires the tail start parameter and a reference to an array of three doubles. 

Modifies the array of three doubles. The first entry is the ratios of all 
particles. The second is the neutron ratio an the third is the photon. 
*/
void mean_qval(string dataFileName, int tailStart, double* ratioArray,
			   TH1D* bothHist, TH1D* neutronHist, TH1D* photonHist);

int main(int argc, char **argv){	
	
	string dataFileName = argv[1];

	double mean_array[3] = {0,0,0};

	TCanvas * canvas = new TCanvas("Tail vs Total Integration", "Both Particles' Tail/Total");
	
	TH1D * bothHist = new TH1D("a_h", "All Particles' Total/Tail", 100, 0.0, 1.5);
	TH1D * neutronHist = new TH1D("n_h", "Neutron Total/Tail", 100, 0.0, 1.5);
	TH1D * photonHist = new TH1D("p_h", "Photon Total/Tail", 100, 0.0, 1.5);

	mean_qval(argv[1], 20, mean_array,
			  bothHist, neutronHist, photonHist);
			  
	neutronHist->Draw();
	
	string plotName = DATA_FOLDER_PATH  + dataFileName + "_qval.pdf";
	canvas->Print((plotName+"[").c_str(),"pdf");
	
}

void mean_qval(string dataFileName, int tailStart,
			   TH1D* bothHist, TH1D* neutronHist, TH1D* photonHist){
	string dataFolderPrefix = "../data/";
	
	TFile* par_file = new TFile((dataFolderPrefix + dataFileName + "_coarse_time.root").c_str(), "UPDATE");
	TTree* par_tree =( TTree*)par_file->Get((dataFileName + "_coarse_time").c_str());
	if (!par_tree){
		cout << "Could not find the " <<  dataFileName << "_coarse_time.root file" << endl;
	}
	
	TFile* raw_file = new TFile((dataFolderPrefix + dataFileName + "_raw_data.root").c_str(), "READ");
	TTree* raw_tree = (TTree*)raw_file->Get(Form("%s_raw", dataFileName.c_str()));
	if (!raw_tree){
			cout << "Could not find the " <<  dataFileName << "raw_time.root file" << endl;
		}

	
	double ped_n, ped_f, amp_n, amp_f, t_max_n, t_max_f, t_50_n, t_50_f;
	
	int is_good_n = -1, 
		is_good_f = -1,
		adc_n[N_width*32] = {0},
		adc_f[N_width*32] = {0};
	
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

	TH1D* qval=new TH1D("qval","Distribution of tail Q/total Q; tail Q/total Q; counts", 100, 0., 2);
	TH1D* qval_n=new TH1D("qval_n","Distribution of tail Q/total Q; tail Q/total Q; counts", 100, 0., 2);
	TH1D* qval_p=new TH1D("qval_p","Distribution of tail Q/total Q; tail Q/total Q; counts", 100, 0., 2);

	int n_neutron=0;
	int n_photon=0;

	cout << "Starting loop";
	for (int i=0; i<n_tot; i++){
		par_tree->GetEntry(i);
		raw_tree->GetEntry(i);

		if (is_good_f  && is_good_n && t_max_f!=0. && t_max_n!=0.){
			double start_time = t_50_f - 20.;
			
			int start_TDC = (int)ceil(start_time);
			if (start_TDC%2){
				start_TDC += 1;
			}
			start_TDC/=2;

			int tail_TDC = start_TDC+20;
			int end_TDC = start_TDC+110;

			double total_area = 0.;
			double tail_area = 0.;

			for (int j = start_TDC; j <tail_TDC; j++){
				total_area += adc_f[j]-ped_f;
			}

			// The integration of the waveforms.
			for (int j=tail_TDC; j<end_TDC; j++){
				total_area += (double)adc_f[j]-ped_f;
				tail_area  += (double)adc_f[j]-ped_f;
			}
			
			double ratio=tail_area/total_area;
			qval->Fill(ratio);
			
			double t_50_diff = t_50_f - t_50_n;
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
	
	bothHist = qval;
	neutronHist = qval_n;
	photonHist = qval_p;

	par_file->Close();
	raw_file->Close();

}

