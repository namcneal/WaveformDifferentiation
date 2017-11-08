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

int channel_norm_fit(int* adc_, double* ped, double* amp, double* t_max, double* t_50){
	double ped_lv=0.;
	for (int i=0; i<10; i++) ped_lv+=(double)adc_[i];
	ped_lv/=10.;
	int TDC_start_peak=0;
	int TDC_end_peak=0;
	int TDC_max=0;
	double ADC_max=0;
	for (int i=0; i<300; i++){
		if ((double)adc_[i]>(ADC_max)) {
			TDC_max=i;
			ADC_max=(double)adc_[i];
		}
		if (((double)adc_[i]>(ped_lv+50.))&&(TDC_start_peak==0)) TDC_start_peak=i-5;
	}
	TDC_end_peak=TDC_max+2;
	if (TDC_start_peak<0) TDC_start_peak=0;
	cout << t_max << endl;
	*t_max = (double)TDC_max * PERTDC; 
	int n_fit=TDC_end_peak-TDC_start_peak+1;
	
//	// /*
	int saturation_flag=0;
	int fardoublepeak_flag=0;
	int closedoublepeak_flag=0;
	int count=0;
	for (int i=TDC_max-3; i<TDC_max+4; i++){
		if (adc_[i]==ADC_max) count++;
	}
	if (count>2) saturation_flag=1;
	for (int i=0; i<TDC_max-10; i++){
		if (adc_[i]>(0.4*(ADC_max-ped_lv)+ped_lv)) {fardoublepeak_flag=1; break;}
	}
	for (int i=TDC_max+30; i<N_width*32; i++){
		if (adc_[i]>(0.4*(ADC_max-ped_lv)+ped_lv)) {fardoublepeak_flag=1; break;}
	}
	for (int i=TDC_max-3; i<TDC_max; i++){
		if (adc_[i]<adc_[i-1]) {closedoublepeak_flag=1; break;}
	}
	for (int i=TDC_max+1; i<TDC_max+3; i++){
		if (adc_[i]<adc_[i+1]) {closedoublepeak_flag=1; break;}
	}
	// */
	
	//if ((n_fit>8)&&(n_fit<15)&&(ped_lv<400.)){
	if ((n_fit>8)&&(n_fit<15)&&(ped_lv<400.)&&(saturation_flag==0)&&(fardoublepeak_flag==0)&&(closedoublepeak_flag==0)){
		double ADC_fit[n_fit];
		double t_fit[n_fit];
		double dy[n_fit];
		for (int i=TDC_start_peak; i<(TDC_end_peak+1); i++){
			ADC_fit[i-TDC_start_peak]=(double)adc_[i];
			t_fit[i-TDC_start_peak]=PERTDC*(double)i;
			//dy[i-TDC_start_peak]=TMath::Sqrt(fabs(((double)adc_[i]-ped_lv)));
			dy[i-TDC_start_peak]=1.;//equal weight put to each
		}
		TGraphErrors* g1_fit=new TGraphErrors(n_fit, t_fit, ADC_fit, 0, dy);
		TF1 *func1=new TF1("f", "[3]+[0]*TMath::Gaus(x, [1], [2], 0)", ((double)TDC_start_peak-0.5)*PERTDC, ((double)TDC_end_peak+0.5)*PERTDC);
		func1->SetParameter(0, ADC_max-ped_lv);
		func1->SetParameter(1, *t_max);
		func1->SetParameter(2, 5.);
		func1->FixParameter(3, ped_lv);
		//func1->SetParameter(3, ped_lv);
		func1->SetParameter(3, ped_lv);
		func1->SetParLimits(0, 0.1*(ADC_max-ped_lv), 10.*(ADC_max-ped_lv));
		func1->SetParLimits(1, *t_max-5.*PERTDC, *t_max+5.*PERTDC);
		func1->SetParLimits(2, 0., 200.);
		//func1->SetParLimits(3, ped_lv-250., ped_lv+250.);
		TFitResultPtr frp1=g1_fit->Fit(func1, "QRSM");
		*ped=func1->GetParameter(3);
		*amp=func1->GetParameter(0);
		*t_50=func1->GetX(func1->GetParameter(3)+0.5*(func1->GetParameter(0)), t_fit[0]-5.*PERTDC, *t_max);
		double chi2=func1->GetChisquare(); 
		//cout<<chi2<<" ";
		//if (chi2<n_fit*3.&&*t_max>200.){
		if (*t_max>200.&&ADC_max<4090.){
			delete g1_fit;
			delete func1;
			return 0;
		}
		else return 2;
	}
	else return 1;
}

int main(int argc, char **argv){
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	string filePathPrefix = argv[1],
		   fileName = argv[2];
	
	
	string new_file_name = filePathPrefix + fileName + "_coarse_time.root",
		   new_data_tree_name = fileName + "_coarse_time";
	TFile* new_file=new TFile(new_file_name.c_str(), "recreate");
	TTree* new_tree=new TTree(new_data_tree_name.c_str(),(fileName+"_coarse_time").c_str());
	if (!new_tree){
		cout << "Could not creat the new " << fileName << "_course_time tree." << endl;
	}
	
	string raw_file_name = filePathPrefix + fileName + "_raw_data.root",
		   raw_data_tree_name = fileName + "_raw";
	TFile* raw_file=new TFile(raw_file_name.c_str(), "READ");
	TTree* raw_tree = (TTree*)raw_file->Get(raw_data_tree_name.c_str());
	if (!raw_tree){
		cout << "\n\nCannot find raw data tree" << endl << endl;
	}
	
	double ped_n = 0, 
		   ped_f = 0, 
		   amp_n = 0, 
		   amp_f = 0, 
		   t_max_n = 0, 
		   t_max_f = 0, 
		   t_50_n = 0, 
		   t_50_f = 0;
	
	int is_good_n = 0, 
		is_good_f = 0;
		
	int adc_n[N_width*32] = {0},
		adc_f[N_width*32] = {0};
	
	cout << endl << "Setting branch addresses for the new tree." << endl;
	new_tree->Branch("ped_n", &ped_n, "ped_n/D");
	new_tree->Branch("ped_f", &ped_f, "ped_f/D");
	new_tree->Branch("amp_n", &amp_n, "amp_n/D");
	new_tree->Branch("amp_f", &amp_f, "amp_f/D");
	new_tree->Branch("t_max_n", &t_max_n, "t_max_n/D");
	new_tree->Branch("t_max_f", &t_max_f, "t_max_f/D");
	new_tree->Branch("t_50_n", &t_50_n, "t_50_n/D");
	new_tree->Branch("t_50_f", &t_50_f, "t_50_f/D");
	
	cout << "Setting branch addresses for the raw tree" << "...adc n";
	raw_tree->SetBranchAddress("adc_n", adc_n);
	cout << "...adc f";
	raw_tree->SetBranchAddress("adc_f", adc_f);
	cout << "...n quality";
	raw_tree->SetBranchAddress("n_quality", &is_good_n);
	cout << "...f quality." << endl;
	raw_tree->SetBranchAddress("f_quality", &is_good_f);
	
	cout << "Geting the total number of entries." << endl;
	int n_tot=raw_tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<fileName<<"_raw_data.root"<<endl;
	
	int n_bad_data=0;
	int n_bad_fit=0;
	int n_fit=0;
	
	TH1D* h_t_diff_50=new TH1D("h_t_diff_50", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);
	TH1D* h_t_diff_max=new TH1D("h_t_diff_max", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);
	
	cout << "Starting loop" << endl;
	for (int i=0; i<n_tot; i++){
		
		cout << "Getting entry from raw data tree;" << endl;
		raw_tree->GetEntry(i);
		
		cout << "Running the normalization function and filling new tree." << endl;
		int flag_n=channel_norm_fit(adc_n, &ped_n, &amp_n, &t_max_n, &t_50_n);
		int flag_f=channel_norm_fit(adc_f, &ped_f, &amp_f, &t_max_f, &t_50_f);
		cout << flag_n << " " << flag_f<<endl;
		new_tree->Fill();
		
		cout << "Determining good and bad data" << endl;
		if (is_good_n==0||is_good_f==0) n_bad_data++;
		else{
			if (flag_n==0&&flag_f==0){
				n_fit++;
				h_t_diff_50->Fill(t_50_f-t_50_n);
				h_t_diff_max->Fill(t_max_f-t_max_n);
			}
			else n_bad_fit++;
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_bad_data<<" events have empty channel(s)."<<endl;
	cout<<"    "<<n_bad_fit<<" events have bad fitting in either channel."<<endl;
	cout<<"    "<<n_fit<<" events have proper fits and time difference of signals were calculated."<<endl;
	new_file->cd();
	new_tree->Write();
	
	string plot_name = filePathPrefix  + fileName + "_coarse_time.pdf";
	TCanvas* c0=new TCanvas("0", "", 800, 600);
	
	gStyle->SetOptStat(110010);
	h_t_diff_50->SetLineColor(kRed);
	h_t_diff_max->SetLineColor(kBlack);
	
	h_t_diff_50->Draw();
	c0->Print((plot_name+"(").c_str(),"pdf");
	
	c0->SetLogy();
	h_t_diff_max->SetMinimum(0.8);
	h_t_diff_max->Draw();
	h_t_diff_50->Draw("SAME");
	TLegend* leg=new TLegend(.6, .9, .9, .7, "");
	leg->SetFillColor(0);
	leg->AddEntry(h_t_diff_50, "50% rising edge");
	leg->AddEntry(h_t_diff_max, "maximum point (discrete)");
	leg->Draw("SAME");
	c0->Print((plot_name+")").c_str(),"pdf");
	
	return 0;
}


