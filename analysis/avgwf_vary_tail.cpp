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

double landaufun(double *x, double *par) {
	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation),
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.

	// Numeric constants
	double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	double mpshift  = -0.22278298;       // Landau maximum location

	// Control constants
	double np = 800.0;      // number of convolution steps
	double sc =   7.0;      // convolution extends to +-sc Gaussian sigmas
	// Variables
	double xx;
	double mpc;
	double fland;
	double sum = 0.0;
	double xlow,xupp;
	double step;
	double i;
	// MP shift correction
	mpc = par[1] - mpshift * par[0];
	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];
	step = (xupp-xlow) / np;
	
	// Convolution integral of Landau and Gaussian by sum
	for(i=1.0; i<=np/2; i++) {
		xx = xlow + (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
	
		xx = xupp - (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
	}
	
	return (par[2] * step * sum * invsq2pi / par[3]);
}

//need to set TF1 and names independently
//need to plot afterwards
void fit_func(TProfile* plot, TF1* func, int n_par, double* par_init, double* par_llim, double* par_hlim, TCanvas* c, int div_no, string const &formula, int use_L=0){
	c->cd(div_no);
	TPad* pad1=new TPad("pad1","",0,0.3,1,1);
	pad1->SetFillColor(0);
	pad1->SetFillStyle(4000);
	pad1->SetFrameFillStyle(0);
	pad1->Draw();
	pad1->cd();
	gStyle->SetOptFit(1111);
	
	plot->SetMarkerColor(kBlue);
	plot->SetLineColor(kBlue);
	func->SetLineColor(kRed);
	
	for (int i=0; i<n_par; i++){
		func->SetParameter(i, par_init[i]);
		func->SetParLimits(i, par_llim[i], par_hlim[i]);
	}
	if (use_L==0) plot->Fit(func,"RQM");
	else plot->Fit(func,"RQLM");
	plot->Draw("E");
	plot->GetXaxis()->SetLimits(-20.,200.);
	plot->GetXaxis()->SetLabelSize(0.06);
	plot->GetYaxis()->SetLabelSize(0.06);
	plot->GetXaxis()->SetTitleSize(0.06);
	plot->GetYaxis()->SetTitleSize(0.06);
	plot->GetXaxis()->SetTitleOffset(0.6);
	plot->GetYaxis()->SetTitleOffset(0.6);
	func->DrawClone("SAME");
	
	TPaveText *pt = new TPaveText(.35,.9-(2.5+n_par)*0.075,.95,.93, "brNDC");
	pt->SetTextSize(0.055);
	pt->AddText(formula.c_str());
	pt->AddText(Form("#chi^{2} / NDF        %.3lf/%i", func->GetChisquare(), func->GetNDF()));
	for (int i=0; i<n_par; i++){
		pt->AddText(Form("%s       %lf #pm %lf", func->GetParName(i), func->GetParameter(i), func->GetParError(i)));
	}
	pt->Draw();
	
	c->cd(div_no);
	TPad* pad2=new TPad("pad2","",0,0,1,0.3);
	pad2->SetFillColor(0);
	pad2->SetFillStyle(4000);
	pad2->SetFrameFillStyle(0);
	pad2->Draw();
	pad2->cd();
	gStyle->SetOptFit(1111);
	
	int n_bin=440;
	double error[n_bin];
	double xline[n_bin];
	double t_fit[n_bin];
	for (int i=1; i<=n_bin; i++){
		t_fit[i]=-20.25+i*0.5;
		error[i]=(plot->GetBinContent(i)-func->Eval(t_fit[i]))/plot->GetBinError(i);
		xline[i]=.25;
	}
	
	TGraphErrors* err=new TGraphErrors(n_bin, t_fit, error, xline, 0);
	err->SetTitle(";;[ADC(t)-f(t)]/#sigma");
	err->SetMarkerColor(kBlue);
	err->SetLineColor(kBlue);
	err->Draw("AP");
	err->GetXaxis()->SetLimits(-20.,200.);
	err->GetXaxis()->SetLabelSize(0.125);
	err->GetYaxis()->SetLabelSize(0.125);
	err->GetXaxis()->SetTitleSize(0.125);
	err->GetYaxis()->SetTitleSize(0.125);
	err->GetYaxis()->SetTitleOffset(0.4);
	TLine* line=new TLine(-20.,0.,200.,0.);
	line->SetLineWidth(2);
	line->SetLineColor(kBlack);
	line->DrawClone("SAME");
	
	return;
}

void calculateQVals(double tailStart, TProfile* avgNeutron, TProfile* avgPhoton,
					double &neutronQVal, double &photonQVal, 
					double &ratioQVal, double &diffQVal);
					  
int main(int argc, char **argv){
	cout << "Reading in the data and creating the trees" << endl;
	string dataFileName = argv[1],
		   dataFolderPath = "../data/";
		   
	string plot_name = dataFolderPath + dataFileName + "PLEASE.pdf";
	
	TFile* par_file=new TFile((dataFolderPath + dataFileName + "_coarse_time.root").c_str(), "UPDATE");
	TTree* par_tree=(TTree*)par_file->Get((dataFileName + "_coarse_time").c_str());
	
	TFile* raw_file=new TFile((dataFolderPath + dataFileName+"_raw_data.root").c_str(), "READ");
	TTree* raw_tree=(TTree*)raw_file->Get((dataFileName + "_raw").c_str());
	
	double ped_n, ped_f, amp_n, amp_f, t_max_n, t_max_f, t_50_n, t_50_f;
	int is_good_n, is_good_f;
	int adc_n[N_width*32];
	int adc_f[N_width*32];
	
	cout << "Setting branch addresses to read into." << endl;
	par_tree->SetBranchAddress("ped_n", &ped_n);
	par_tree->SetBranchAddress("ped_f", &ped_f);
	par_tree->SetBranchAddress("amp_n", &amp_n);
	par_tree->SetBranchAddress("amp_f", &amp_f);
	par_tree->SetBranchAddress("t_max_n", &t_max_n);
	par_tree->SetBranchAddress("t_max_f", &t_max_f);
	par_tree->SetBranchAddress("t_50_n", &t_50_n);
	par_tree->SetBranchAddress("t_50_f", &t_50_f);
	
	raw_tree->SetBranchAddress("adc_n", adc_n);
	raw_tree->SetBranchAddress("adc_f", adc_f);
	raw_tree->SetBranchAddress("n_quality", &is_good_n);
	raw_tree->SetBranchAddress("f_quality", &is_good_f);
	
	int n_tot=raw_tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<dataFileName<<"_raw_data.root"<<endl;
	
//	TGraph* sct_pts_n=new TGraph();
//	sct_pts_n->SetTitle("Normalized Waveform: Neutron; time [ns]; Percentage");
	TH2D* hist_n=new TH2D("hist_n","Normalized Waveform: Neutron; time [ns]; Percentage",440,-20., 200.,140,-.2,1.2);
	
//	TGraph* sct_pts_p=new TGraph();
//	sct_pts_p->SetTitle("Normalized Waveform: Photon; time [ns]; Percentage");
	TH2D* hist_p=new TH2D("hist_p","Normalized Waveform: Photon; time [ns]; Percentage",440,-20., 200.,140,-.2,1.2);
	
	int n_neutron=0;
	int n_photon=0;
	
	for (int i=0; i<n_tot; i++){
		par_tree->GetEntry(i);
		raw_tree->GetEntry(i);
		
		if (is_good_f==1 && is_good_n==1 && t_max_f!=0. && t_max_n!=0.){
			double t_50_diff=t_50_f-t_50_n;
			double start_time=t_50_f-20.;
			int start_TDC=(int)ceil(start_time);
			if (start_TDC%2) start_TDC+=1;
			start_TDC/=2;
			
			if (t_50_diff>neutron_low && t_50_diff<neutron_high){
				for (int j=0; j<110; j++){
//					sct_pts_n->SetPoint(n_neutron*110+j, 2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
					hist_n->Fill(2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
				}
				n_neutron++;
			}
			if (t_50_diff>photon_low && t_50_diff<photon_high){
				for (int j=0; j<110; j++){
//					sct_pts_p->SetPoint(n_photon*110+j, 2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
					hist_p->Fill(2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
				}
				n_photon++;
			}
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_neutron<<" neutron events are present."<<endl;
	cout<<"    "<<n_photon<<" photon events are present."<<endl;
	
	TProfile* avg_n=hist_n->ProfileX();
	TProfile* avg_p=hist_p->ProfileX();
	
	Int_t numSteps = 150,
		  minTailStart = 10,
		  maxTailStart = 180,
		  minTailEnd = 170,
		  maxTailEnd = 170;

	TArrayI params(5);
	params.AddAt(numSteps, 0);
	params.AddAt(minTailStart, 1);
	params.AddAt(maxTailStart, 2);
	params.AddAt(minTailEnd, 3);
	params.AddAt(maxTailEnd, 4);


	double stepSize = (maxTailStart - minTailStart) / numSteps;
		   
	double tailStart,
		   neutronQVal,
		   photonQVal,
		   photonNeutronQValRatio,
		   photonNeutronQValDiff;
		   
	TFile* variedTailFile = new TFile((dataFolderPath + dataFileName + "_varied_tail_qvals.root").c_str(), "RECREATE");
	
	TTree* qValTree = new TTree("data", "Varied Tail Q-value Results");
	qValTree->Branch("tailStart", &tailStart);
	qValTree->Branch("neutronQVal", &neutronQVal);
	qValTree->Branch("photonQVal", &photonQVal);
	qValTree->Branch("ratioQVal", &photonNeutronQValRatio);
	qValTree->Branch("diffQVal", &photonNeutronQValDiff);

	avg_n->Write();
	avg_p->Write();
	variedTailFile->WriteObjectAny(&params,"TArrayI","calculationParams");
	
	for(tailStart = minTailStart; tailStart <= maxTailStart; tailStart += stepSize){
		calculateQVals(tailStart, avg_n, avg_p,
					   neutronQVal, photonQVal, 
					   photonNeutronQValRatio, photonNeutronQValDiff);
		qValTree->Fill();
	}
	qValTree->Write();
	variedTailFile->Close();
	
	
	TFile* readFile=new TFile((dataFolderPath + dataFileName+"_varied_tail_qvals.root").c_str(), "READ");
	TTree* readTree=(TTree*)readFile->Get("data");

	double tailStarts[numSteps], 
		   neutronQVals[numSteps],
		   photonQVals[numSteps],
		   photonNeutronQValRatios[numSteps],
		   photonNeutronQValDiffs[numSteps];
	
	
	for(int i =0; i <= numSteps; i++){
		readTree->SetBranchAddress("tailStart", tailStarts + i);
		readTree->SetBranchAddress("neutronQVal", neutronQVals + i);
		readTree->SetBranchAddress("photonQVal", photonQVals + i);
		readTree->SetBranchAddress("ratioQVal", photonNeutronQValRatios + i);
		readTree->SetBranchAddress("diffQVal", photonNeutronQValDiffs + i);
		
		readTree->GetEntry(i);
	}
	
	TCanvas* c0=new TCanvas("0", "", 600, 500);
	c0->Divide(1,2);
	
	c0->cd(1);
	TGraph* neutronPlot = new TGraph(numSteps, tailStarts, neutronQVals);
	neutronPlot->SetTitle("Photon (Blue) & Neutron (Red) Individual Q-values and Difference");
	neutronPlot->GetXaxis()->SetLimits(0, 200);
	neutronPlot->SetMinimum(0);
	neutronPlot->SetMaximum(.6);
	neutronPlot->SetMarkerStyle(24);
	neutronPlot->SetMarkerSize(.5);
	neutronPlot->SetMarkerColor(2);
	neutronPlot->SetLineColor(2);
	neutronPlot->Draw("Al");
	
	TGraph* photonPlot = new TGraph(numSteps, tailStarts, photonQVals);
	photonPlot->SetMarkerStyle(24);
	photonPlot->SetMarkerColor(4);
	photonPlot->SetMarkerSize(.5);
	photonPlot->SetLineColor(4);
	photonPlot->Draw("l");
	
	TGraph* ratioPlot = new TGraph(numSteps, tailStarts, photonNeutronQValDiffs);
	ratioPlot->SetMarkerStyle(24);
	ratioPlot->SetMarkerColor(3);
	ratioPlot->SetMarkerSize(.5);
	ratioPlot->SetLineColor(3);
	ratioPlot->Draw("l");
	
	TLine* line_top =new TLine(20,0,20,1.1);
	line_top->SetLineColor(kBlack);
	line_top->Draw("SAME");
		
	c0->cd(2);
	avg_n->SetLineColor(kRed);
	avg_p->SetLineColor(kBlue);
	avg_n->SetLineWidth(3);
	avg_p->SetLineWidth(3);
	avg_n->SetStats(false);
	avg_p->SetStats(false);
	avg_n->GetXaxis()->SetLimits(0, 200);
	avg_p->GetXaxis()->SetLimits(0, 200);
	avg_n->Draw("HIST C");
	avg_p->Draw("A HIST SAME C");
	
	TLine* line_bottom =new TLine(20,0,20,1.1);
	line_bottom->SetLineColor(kBlack);
	line_bottom->Draw("SAME");
	
	
	c0->Modified();
	c0->Update();
	c0->SaveAs((dataFolderPath + dataFileName + "_varied_qvals_DIFFERENCE.png").c_str(), "png");
	
	par_file->Close();
	raw_file->Close();
	return 0;
}

void calculateQVals(double tailStart, TProfile* avgNeutron, TProfile* avgPhoton,
					double &neutronQVal, double &photonQVal, double &ratioQVal, double &diffQVal){
						  
	int tail_start_bin=avgNeutron->FindBin(tailStart);
	double a_main_n=0.,
		   a_tail_n=0.,
		   a_main_p=0.,
		   a_tail_p=0.;
	
	for (int i=1; i<tail_start_bin; i++){
		a_main_n+=avgNeutron->GetBinContent(i);
		a_main_p+=avgPhoton->GetBinContent(i);
	}
	for (int i=tail_start_bin; i<=avgNeutron->FindBin(199.75); i++){
		a_tail_n+=avgNeutron->GetBinContent(i);
		a_tail_p+=avgPhoton->GetBinContent(i);
	}
	
	double neutronCurrentQVal = a_tail_n/(a_tail_n+a_main_n),
		   photonCurrentQVal = a_tail_p/(a_tail_p+a_main_p);
	neutronQVal = neutronCurrentQVal;
	photonQVal  = photonCurrentQVal;
	ratioQVal   = photonCurrentQVal / neutronCurrentQVal;
	diffQVal = abs(photonCurrentQVal - neutronCurrentQVal);
							  
}