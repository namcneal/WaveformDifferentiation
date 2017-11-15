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

int main(int argc, char **argv){
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	const Int_t NRGBs = 5;
	const Int_t NCont = 256;
	Double_t stops[NRGBs] = { 0.00, 0.30, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs] = { 0.00, 0.00, 0.57, 0.90, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.65, 0.95, 0.20, 0.00 };
	Double_t blue[NRGBs] = { 0.51, 0.55, 0.15, 0.00, 0.10 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont); 
	
	cout << "Reading in the data and creating the trees" << endl;
	string dataFileName = argv[1],
		   dataFolderPath = "../data/";
		   
	string plot_name = dataFolderPath + dataFileName + "_avg_wf.pdf";
	
	TFile* par_file=new TFile((dataFolderPath + dataFileName + "_coarse_time.root").c_str(), "UPDATE");
	TTree* par_tree=(TTree*)par_file->Get((dataFileName + "_coarse_time").c_str());
	
	TFile* raw_file=new TFile((dataFolderPath + dataFileName+"_raw_data.root").c_str(), "READ");
	TTree* raw_tree=(TTree*)raw_file->Get((dataFileName + "_raw").c_str());
	
	TCanvas* c0=new TCanvas("0", "", 3200, 1200);
	c0->Divide(2,1);
	c0->Print((plot_name+"[").c_str(),"pdf");
	
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
	
	TGraph* sct_pts_n=new TGraph();
	sct_pts_n->SetTitle("Normalized Waveform: Neutron; time [ns]; Percentage");
	TH2D* hist_n=new TH2D("hist_n","Normalized Waveform: Neutron; time [ns]; Percentage",440,-20., 200.,140,-.2,1.2);
	
	TGraph* sct_pts_p=new TGraph();
	sct_pts_p->SetTitle("Normalized Waveform: Photon; time [ns]; Percentage");
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
					sct_pts_n->SetPoint(n_neutron*110+j, 2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
					hist_n->Fill(2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
				}
				n_neutron++;
			}
			if (t_50_diff>photon_low && t_50_diff<photon_high){
				for (int j=0; j<110; j++){
					sct_pts_p->SetPoint(n_photon*110+j, 2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
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
	avg_n->SetTitle("X-Profile of Normalized Waveform: Neutron; time [ns]; Percentage");
	TProfile* avg_p=hist_p->ProfileX();
	avg_p->SetTitle("X-Profile of Normalized Waveform: Photon; time [ns]; Percentage");
	
	sct_pts_n->SetMarkerStyle(6);
	sct_pts_n->SetMarkerColor(kRed);
	sct_pts_p->SetMarkerStyle(6);
	sct_pts_p->SetMarkerColor(kRed);
	c0->cd(1);
	sct_pts_n->Draw("AP");
	c0->cd(2);
	sct_pts_p->Draw("AP");
	//c0->Print((plot_name).c_str(),"pdf");
	
	c0->cd(1);
	hist_n->Draw("COLZ1");
	c0->cd(2);
	hist_p->Draw("COLZ1");
	c0->Print((plot_name).c_str(),"pdf");
	
	avg_n->SetMarkerColor(kRed);
	avg_p->SetMarkerColor(kBlue);
	avg_n->SetLineColor(kRed);
	avg_p->SetLineColor(kBlue);
	c0->cd(1);
	avg_n->Draw("E");
	c0->cd(2);
	avg_p->Draw("E");
	c0->Print((plot_name).c_str(),"pdf");
	cout<<"Average waveform constructed..."<<endl;
	
	TCanvas* c1=new TCanvas("1", "", 3200, 1200);
	c1->Divide(2,1);
	c1->cd(1);
	avg_n->SetTitle("X-Profile of Normalized Waveform; time [ns]; Percentage");
	avg_p->SetTitle("X-Profile of Normalized Waveform; time [ns]; Percentage");
	avg_n->Draw();
	avg_p->Draw("SAME");
	TLine* line=new TLine(tail_start,0.,tail_start,1.);
	line->SetLineColor(kBlack);
	line->Draw("SAME");
	TLegend* leg=new TLegend(.6, .9, .9, .7, "");
	leg->SetFillColor(0);
	leg->AddEntry(avg_n, "neutron");
	leg->AddEntry(avg_p, "photon");
	leg->AddEntry(line, Form("t=%.1lf ns", tail_start));
	leg->DrawClone("SAME");
	avg_n->SetTitle("X-Profile of Normalized Waveform: Neutron; time [ns]; Percentage");
	avg_p->SetTitle("X-Profile of Normalized Waveform: Photon; time [ns]; Percentage");
	
	c1->cd(2);
	int tail_start_bin=avg_n->FindBin(tail_start);
	double a_main_n=0.;
	double a_tail_n=0.;
	double a_main_p=0.;
	double a_tail_p=0.;
	for (int i=1; i<tail_start_bin; i++){
		a_main_n+=avg_n->GetBinContent(i);
		a_main_p+=avg_p->GetBinContent(i);
	}
	for (int i=tail_start_bin; i<=avg_n->FindBin(199.75); i++){
		a_tail_n+=avg_n->GetBinContent(i);
		a_tail_p+=avg_p->GetBinContent(i);
	}
	TPaveText *pt = new TPaveText(.05,.1,.95,.8);
	pt->AddText(Form("Tail is defined as the part of the waveform to the right of t=%.1lf ns", tail_start));
	pt->AddText("Defining the Tail Area Ratio Q");
	pt->AddText("Q = #frac{Area under the tail}{Total area under the waveform}");
	pt->AddLine(.0,.4,1.,.4);
	pt->AddText(Form("For the averaged photon waveform, Q=%lf", a_tail_p/(a_tail_p+a_main_p)));
	pt->AddText(Form("For the averaged neutron waveform, Q=%lf", a_tail_n/(a_tail_n+a_main_n)));
	pt->Draw();
	
	c1->Print((plot_name).c_str(),"pdf");
	cout<<"Q-value method applied..."<<endl;
	
	//fit func
	//1. simple land.
	TCanvas* c2=new TCanvas("2", "", 3200, 1200);
	gStyle->SetOptFit(1111);
	c2->Divide(2,1);
	
	string form1="A*NormalizedLandau(x, #mu, #sigma)";
	
	TF1* simp_land_n=new TF1("simp_land_n", "[0]*TMath::Landau(x,[1],[2],1)",fit_left,fit_right);
	simp_land_n->SetParNames("A","#mu","#sigma");
	double par_init_1_n[3]={1.,0.,3.};
	double par_llim_1_n[3]={0.,-20.,0.0001};
	double par_hlim_1_n[3]={100000.,200.,30.};
	fit_func(avg_n, simp_land_n, 3, par_init_1_n, par_llim_1_n, par_hlim_1_n, c2, 1, form1);
	
	TF1* simp_land_p=new TF1("simp_land_p", "[0]*TMath::Landau(x,[1],[2],1)",fit_left,fit_right);
	simp_land_p->SetParNames("A","#mu","#sigma");
	double par_init_1_p[3]={1.,0.,3.};
	double par_llim_1_p[3]={0.,-20.,0.0001};
	double par_hlim_1_p[3]={100000.,200.,30.};
	fit_func(avg_p, simp_land_p, 3, par_init_1_p, par_llim_1_p, par_hlim_1_p, c2, 2, form1);
	
	c2->Print((plot_name).c_str(),"pdf");
	cout<<"Simple Landau fit completed..."<<endl;
	
	//fit func
	//2. gaus. convoluted land.
	TCanvas* c3=new TCanvas("3", "", 3200, 1200);
	gStyle->SetOptFit(1111);
	c3->Divide(2,1);
	
	string form2="A*NGaus(#sigma_{gaus})[conv]NLandau(x, #mu, #sigma_{land})";
	
	TF1 *landaus_n=new TF1("landaus_n",landaufun,fit_left,fit_right,4);
	landaus_n->SetParNames("#sigma_{land}","#mu","A","#sigma_{gaus}");
	double par_init_2_n[4]={0.1,0.,100, 0.1};
	double par_llim_2_n[4]={0.0001,-20.,0.,0.0001};
	double par_hlim_2_n[4]={30.,200.,100000.,30.};
	fit_func(avg_n, landaus_n, 4, par_init_2_n, par_llim_2_n, par_hlim_2_n, c3, 1, form2);
	
	TF1 *landaus_p=new TF1("landaus_p",landaufun,fit_left,fit_right,4);
	landaus_p->SetParNames("#sigma_{land}","#mu","A","#sigma_{gaus}");
	double par_init_2_p[4]={0.1,0.,100, 0.1};
	double par_llim_2_p[4]={0.0001,-20.,0.,0.0001};
	double par_hlim_2_p[4]={30.,200.,100000.,30.};
	fit_func(avg_p, landaus_p, 4, par_init_2_p, par_llim_2_p, par_hlim_2_p, c3, 2, form2);
	
	c3->Print((plot_name).c_str(),"pdf");
	cout<<"Gaussian convoluted Landau fit completed..."<<endl;
	
	//fit func
	//3. exp. modified gaus. (one from group meeting)
	TCanvas* c4=new TCanvas("4", "", 3200, 1200);
	gStyle->SetOptFit(1111);
	c4->Divide(2,1);
	
	//string form3="A #frac{#lambda}{2} Exp[#frac{#lambda}{2}(2#mu +#lambda#sigma^{2}-2x)]Erfc[#frac{#mu +#lambda#sigma^{2}-x}{#sqrt{2}#sigma}]";
	string form3="A*Exp[#frac{1}{2}(#frac{#sigma}{#tau})^{2} - #frac{x-#mu}{#tau}]Erfc[#frac{1}{#sqrt{2}}(#frac{#sigma}{#tau} - #frac{x-#mu}{#sigma})]";
	/*
	TF1* emg_n=new TF1("emg_n", "[0]*[1]/2.*TMath::Exp([1]/2.*(2.*[2]+[1]*[3]*[3]-2.*x))*TMath::Erfc(([2]+[1]*[3]*[3]-x)/(TMath::Sqrt(2)*[3]))",fit_left,fit_right);
	emg_n->SetParNames("A","#lambda","#mu","#sigma");
	double par_init_3_n[4]={100.,0.5,0.,1.};
	double par_llim_3_n[4]={0.,0.,-20.,1.};
	double par_hlim_3_n[4]={100000.,10.,200.,30.};
	*/
	TF1* emg_n=new TF1("emg_n", "[0]*TMath::Exp(.5*[3]*[3]/[1]/[1]-(x-[2])/[1])*TMath::Erfc(1/TMath::Sqrt(2)*([3]/[1]-(x-[2])/[3]))",fit_left,fit_right);
	emg_n->SetParNames("A","#tau","#mu","#sigma");
	double par_init_3_n[4]={100.,100.,0.,1.};
	double par_llim_3_n[4]={0.,5.,-20.,0.};
	double par_hlim_3_n[4]={100000.,100.,200.,30.};
	
	fit_func(avg_n, emg_n, 4, par_init_3_n, par_llim_3_n, par_hlim_3_n, c4, 1, form3);
	/*
	TF1* emg_p=new TF1("emg_p", "[0]*[1]/2.*TMath::Exp([1]/2.*(2.*[2]+[1]*[3]*[3]-2.*x))*TMath::Erfc(([2]+[1]*[3]*[3]-x)/(TMath::Sqrt(2)*[3]))",fit_left,fit_right);
	emg_p->SetParNames("A","#lambda","#mu","#sigma");
	double par_init_3_p[4]={100.,0.5,0.,1.};
	double par_llim_3_p[4]={0.,0.,-20.,1.};
	double par_hlim_3_p[4]={100000.,10.,200.,30.};
	*/
	TF1* emg_p=new TF1("emg_p", "[0]*TMath::Exp(.5*[3]*[3]/[1]/[1]-(x-[2])/[1])*TMath::Erfc(1/TMath::Sqrt(2)*([3]/[1]-(x-[2])/[3]))",fit_left,fit_right);
	emg_p->SetParNames("A","#tau","#mu","#sigma");
	double par_init_3_p[4]={100.,100.,0.,1.};
	double par_llim_3_p[4]={0.,5.,-20.,0.};
	double par_hlim_3_p[4]={100000.,100.,200.,30.};
	
	fit_func(avg_p, emg_p, 4, par_init_3_p, par_llim_3_p, par_hlim_3_p, c4, 2, form3);
	
	c4->Print((plot_name).c_str(),"pdf");
	cout<<"Exponentially modified Gaussian fit completed..."<<endl;
	
	
	//fit func
	//4. gumbel
	TCanvas* c5=new TCanvas("5", "", 3200, 1200);
	gStyle->SetOptFit(1111);
	c5->Divide(2,1);
	
	string form4="A Exp[-(x-#mu)/#beta-Exp(-(x-#mu)/#beta)]";
	
	TF1* gumb_n=new TF1("gumb_n","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))",fit_left,fit_right);
	gumb_n->SetParNames("A","#mu","#beta");
	double par_init_4_n[3]={10.,0.,1.};
	double par_llim_4_n[3]={0.,-20.,1.};
	double par_hlim_4_n[3]={100.,200.,30.};
	
	fit_func(avg_n, gumb_n, 3, par_init_4_n, par_llim_4_n, par_hlim_4_n, c5, 1, form4);
	
	TF1* gumb_p=new TF1("gumb_p","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))",fit_left,fit_right);
	gumb_p->SetParNames("A","#mu","#beta");
	double par_init_4_p[3]={10.,0.,1.};
	double par_llim_4_p[3]={0.,-20.,1.};
	double par_hlim_4_p[3]={100.,200.,30.};
	fit_func(avg_p, gumb_p, 3, par_init_4_p, par_llim_4_p, par_hlim_4_p, c5, 2, form4);
	
	c5->Print((plot_name).c_str(),"pdf");
	cout<<"Gumbel fit completed..."<<endl;
	
	//fit func
	//5. double EMG
	TCanvas* c6=new TCanvas("6", "", 3200, 1200);
	gStyle->SetOptFit(1111);
	c6->Divide(2,1);
	
	string form5="r*EMG(A, #tau_{1}, #mu, #sigma)+(1-r)*EMG(A, #tau_{2}, #mu, #sigma)";
	
	TF1* demg_n=new TF1("demg_n", "[0]/[2]*[1]*TMath::Exp(.5*[5]*[5]/[2]/[2]-(x-[4])/[2])*TMath::Erfc(1./TMath::Sqrt(2.)*([5]/[2]-(x-[4])/[5]))+[0]/[3]*(1.-[1])*TMath::Exp(.5*[5]*[5]/[3]/[3]-(x-[4])/[3])*TMath::Erfc(1./TMath::Sqrt(2.)*([5]/[3]-(x-[4])/[5]))",fit_left,fit_right_long);
	demg_n->SetParNames("A","r","#tau_{1}","#tau_{2}","#mu","#sigma");
	double par_init_5_n[6]={1.,.5,10.,40.,0.,1.};
	double par_llim_5_n[6]={0.,0.,5.,20.,-10.,0.};
	double par_hlim_5_n[6]={20.,1.,20.,100.,20.,30.};
	fit_func(avg_n, demg_n, 6, par_init_5_n, par_llim_5_n, par_hlim_5_n, c6, 1, form5);
	
	TF1* demg_p=new TF1("demg_p", "[0]/[2]*[1]*TMath::Exp(.5*[5]*[5]/[2]/[2]-(x-[4])/[2])*TMath::Erfc(1./TMath::Sqrt(2.)*([5]/[2]-(x-[4])/[5]))+[0]/[3]*(1.-[1])*TMath::Exp(.5*[5]*[5]/[3]/[3]-(x-[4])/[3])*TMath::Erfc(1./TMath::Sqrt(2.)*([5]/[3]-(x-[4])/[5]))",fit_left,fit_right_long);
	demg_p->SetParNames("A","r","#tau_{1}","#tau_{2}","#mu","#sigma");
	double par_init_5_p[6]={1.,.5,10.,40.,0.,1.};
	double par_llim_5_p[6]={0.,0.,5.,20.,-10.,0.};
	double par_hlim_5_p[6]={20.,1.,20.,100.,20.,30.};
	fit_func(avg_p, demg_p, 6, par_init_5_p, par_llim_5_p, par_hlim_5_p, c6, 2, form5);
	
	c6->Print((plot_name).c_str(),"pdf");
	cout<<"Double exponentially modified Gaussian fit completed..."<<endl;
	
	par_file->cd();
	if (par_file->GetListOfKeys()->Contains("hist_n")) par_file->Delete("hist_n");
	if (par_file->GetListOfKeys()->Contains("hist_p")) par_file->Delete("hist_p");
	hist_n->Write();
	hist_p->Write();
	
	c0->Print((plot_name+"]").c_str(),"pdf");
	par_file->Close();
	raw_file->Close();
	return 0;
}
