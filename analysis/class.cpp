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
# include <map>
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
# include "TColor.h"

using namespace std;
const double CHI2THRESHOLD_N=100;
const double CHI2THRESHOLD_F=100;
const double PERTDC=2.;
const int N_width=15;
const double singlePE=14.;

///////////////////////////////////////////
// use as "analysis <filename>"
///////////////////////////////////////////

void bichannel_plot(string const &pdfout, double* adc_n, double* adc_f, double* time, int sampling_len, int index){
	TCanvas* c=new TCanvas("", "", 800, 600);
	TGraph* g1=new TGraph(sampling_len, time, adc_n);
	TGraph* g2=new TGraph(sampling_len, time, adc_f);
	TPad *pad1 = new TPad("pad1","",0,0,1,1);
	TPad *pad2 = new TPad("pad2","",0,0,1,1);
	pad2->SetFillColor(0);
	pad2->SetFillStyle(4000);
	pad2->SetFrameFillStyle(0);
	pad1->Draw();
	pad1->cd();
	g1->SetMarkerStyle(24);
	g1->SetMarkerColor(kRed);
	g1->SetLineColor(kRed);
	g2->SetMarkerStyle(25);
	g2->SetMarkerColor(kBlue);
	g2->SetLineColor(kBlue);
	g1->SetTitle("; Time[ns]; ADC(Red)");
	g2->SetTitle(Form("Instance #%i; Time[ns]; ADC(Blue)", index+1));
	g1->SetMarkerSize(.7);
	g1->SetLineWidth(2);
	g1->SetFillColor(0);
	g2->SetMarkerSize(.7);
	g2->SetLineWidth(2);
	g2->SetFillColor(0);
	g1->Draw("ALP");
	g1->GetYaxis()->SetTitleOffset(1.4);
	pad2->Draw();
	pad2->cd();
	g2->Draw("ALPY+");
	g2->GetYaxis()->SetTitleOffset(1.4);
	TLegend* leg=new TLegend(.5, .8, .8, .6, "");
	leg->SetFillColor(0);
	leg->AddEntry(g1, "near site");
	leg->AddEntry(g2, "far site");
	leg->Draw("SAME");
	c->Print(pdfout.c_str(), "pdf");
	delete c;
}

class channel{
private:
	int sampling_len_;
	int adc_[N_width*32];
	double t_max_;
	double ADC_max_;
	double t_50_;
	double t_50_ref_;
	
	double fp_A_0_;
	
	double fp_demg_A_;
	double fp_demg_r_;
	double fp_demg_tau1_;
	double fp_demg_tau2_;
	double fp_demg_mu_;
	double fp_demg_sigma_;
	double fp_demg_chi2_;
	int is_good_fit_demg_;
	
	double fp_land_A_;
	double fp_land_mu_;
	double fp_land_sigma_;
	double fp_land_chi2_;
	int is_good_fit_land_;
	
public:
	channel()=default;
	void init(){
		sampling_len_=N_width*32;
		for (int i=0; i<sampling_len_; i++) adc_[i]=0;
		t_max_=0.;
		ADC_max_=0.;
		t_50_=0.;
		t_50_ref_=0.;
		
		fp_A_0_=-999.;
		
		fp_demg_A_=0.;
		fp_demg_r_=0.;
		fp_demg_tau1_=0.;
		fp_demg_tau2_=0.;
		fp_demg_mu_=0.;
		fp_demg_sigma_=0.;
		fp_demg_chi2_=0.;
		is_good_fit_demg_=0;
		
		fp_land_A_=0.;
		fp_land_mu_=0.;
		fp_land_sigma_=0.;
		fp_land_chi2_=0.;
		is_good_fit_land_=0;
		
		return;
	}
	int sampling_len(){return sampling_len_;}
	int* adc(){return adc_;}
	int adc_i(int i){return adc_[i];}
	double t_max(){return t_max_;}
	double ADC_max(){return ADC_max_;}
	double t_50(){return t_50_;}
	double t_50_ref(){return t_50_ref_;}
	
	double fp_A_0(){return fp_A_0_;}
	
	double fp_demg_A(){return fp_demg_A_;}
	double fp_demg_r(){return fp_demg_r_;}
	double fp_demg_tau1(){return fp_demg_tau1_;}
	double fp_demg_tau2(){return fp_demg_tau2_;}
	double fp_demg_mu(){return fp_demg_mu_;}
	double fp_demg_sigma(){return fp_demg_sigma_;}
	double fp_demg_chi2(){return fp_demg_chi2_;}
	int is_good_fit_demg(){return is_good_fit_demg_;}
	double fp_land_A(){return fp_land_A_;}
	double fp_land_mu(){return fp_land_mu_;}
	double fp_land_sigma(){return fp_land_sigma_;}
	double fp_land_chi2(){return fp_land_chi2_;}
	int is_good_fit_land(){return is_good_fit_land_;}
	
	void setProcessedBranches(TTree* processed, char label){//label=='n'/'f'
		processed->Branch(Form("t_max_%c", label), &t_max_, Form("t_max_%c/D", label));
		processed->Branch(Form("ADC_max_%c", label), &ADC_max_, Form("ADC_max_%c/D", label));
		processed->Branch(Form("t_50_%c", label), &t_50_, Form("t_50_%c/D", label));
		processed->Branch(Form("fp_A_0_%c", label), &fp_A_0_, Form("fp_A_0_%c/D", label));
		processed->Branch(Form("fp_demg_A_%c", label), &fp_demg_A_, Form("fp_demg_A_%c/D", label));
		processed->Branch(Form("fp_demg_r_%c", label), &fp_demg_r_, Form("fp_demg_r_%c/D", label));
		processed->Branch(Form("fp_demg_tau1_%c", label), &fp_demg_tau1_, Form("fp_demg_tau1_%c/D", label));
		processed->Branch(Form("fp_demg_tau2_%c", label), &fp_demg_tau2_, Form("fp_demg_tau2_%c/D", label));
		processed->Branch(Form("fp_demg_mu_%c", label), &fp_demg_mu_, Form("fp_demg_mu_%c/D", label));
		processed->Branch(Form("fp_demg_sigma_%c", label), &fp_demg_sigma_, Form("fp_demg_sigma_%c/D", label));
		processed->Branch(Form("fp_demg_chi2_%c", label), &fp_demg_chi2_, Form("fp_demg_chi2_%c/D", label));
		processed->Branch(Form("is_good_fit_demg_%c", label), &is_good_fit_demg_, Form("is_good_fit_demg_%c/I", label));
		processed->Branch(Form("fp_land_A_%c", label), &fp_land_A_, Form("fp_land_A_%c/D", label));
		processed->Branch(Form("fp_land_mu_%c", label), &fp_land_mu_, Form("fp_land_mu_%c/D", label));
		processed->Branch(Form("fp_land_sigma_%c", label), &fp_land_sigma_, Form("fp_land_sigma_%c/D", label));
		processed->Branch(Form("fp_land_chi2_%c", label), &fp_land_chi2_, Form("fp_land_chi2_%c/D", label));
		processed->Branch(Form("is_good_fit_land_%c", label), &is_good_fit_land_, Form("is_good_fit_land_%c/I", label));
		return;
	}
	void setRawBranchAddresses(TTree* raw, char label){//label=='n'/'f'
		raw->SetBranchAddress(Form("adc_%c", label), adc_);
		return;
	}
	void setCoarseBranchAddresses(TTree* coarse, char label){//label=='n'/'f'
		coarse->SetBranchAddress(Form("t_50_%c", label), &t_50_ref_);
		return;
	}
	void setProcessedBranchAddresses(TTree* processed, char label){//label=='n'/'f'
		processed->SetBranchAddress(Form("t_max_%c", label), &t_max_);
		processed->SetBranchAddress(Form("ADC_max_%c", label), &ADC_max_);
		processed->SetBranchAddress(Form("t_50_%c", label), &t_50_);
		processed->SetBranchAddress(Form("fp_A_0_%c", label), &fp_A_0_);
		processed->SetBranchAddress(Form("fp_demg_A_%c", label), &fp_demg_A_);
		processed->SetBranchAddress(Form("fp_demg_r_%c", label), &fp_demg_r_);
		processed->SetBranchAddress(Form("fp_demg_tau1_%c", label), &fp_demg_tau1_);
		processed->SetBranchAddress(Form("fp_demg_tau2_%c", label), &fp_demg_tau2_);
		processed->SetBranchAddress(Form("fp_demg_mu_%c", label), &fp_demg_mu_);
		processed->SetBranchAddress(Form("fp_demg_sigma_%c", label), &fp_demg_sigma_);
		processed->SetBranchAddress(Form("fp_demg_chi2_%c", label), &fp_demg_chi2_);
		processed->SetBranchAddress(Form("is_good_fit_demg_%c", label), &is_good_fit_demg_);
		processed->SetBranchAddress(Form("fp_land_A_%c", label), &fp_land_A_);
		processed->SetBranchAddress(Form("fp_land_mu_%c", label), &fp_land_mu_);
		processed->SetBranchAddress(Form("fp_land_sigma_%c", label), &fp_land_sigma_);
		processed->SetBranchAddress(Form("fp_land_chi2_%c", label), &fp_land_chi2_);
		processed->SetBranchAddress(Form("is_good_fit_land_%c", label), &is_good_fit_land_);
		return;
	}
	void fit_func(int func_mode, int index, char label, string const &pdfout, int plotflag=0){//plotflag:2: by chi2(bad fit)*****func_mode: 1=demg; 2=landau
		//extract N_max, ADC_max, find fitting range
		string labelname="";
		double CHI2THRESHOLD=0.;
		if (label=='f') {labelname+="far"; CHI2THRESHOLD=CHI2THRESHOLD_F;} 
		if (label=='n') {labelname+="near"; CHI2THRESHOLD=CHI2THRESHOLD_N;}
		double ped_lv=0.;
		for (int i=0; i<10; i++) ped_lv+=(double)adc_[i];
		ped_lv/=10.;
		 /* check and fix random large point
		for (int i=1; i<sampling_len_-1; i++){
			if (adc_[i-1]<(ped_lv+20.)&&adc_[i+1]<(ped_lv+20.)&&adc_[i]>500.) adc_[i]=(adc_[i-1]+adc_[i+1])/2;
		}
		 */
		int TDC_start_peak=0;
		int TDC_end_peak=0;
		int TDC_max=0;
		ADC_max_=0.;
		for (int i=0; i<300; i++){
			if ((double)adc_[i]>(ADC_max_)) {
				TDC_max=i;
				ADC_max_=(double)adc_[i];
			}
			if (((double)adc_[i]>(ped_lv+50.))&&(TDC_start_peak==0)) TDC_start_peak=i-11;
		}
		TDC_end_peak=TDC_max+75;
		if (TDC_start_peak<0) TDC_start_peak=0;
		t_max_=(double)TDC_max*PERTDC;
		int n_fit=TDC_end_peak-TDC_start_peak+1;
		
		int saturation_flag=0;
		int fardoublepeak_flag=0;
		int closedoublepeak_flag=0;
		int count=0;
		for (int i=TDC_max-3; i<TDC_max+4; i++){
			if (adc_[i]==ADC_max_) count++;
		}
		if (count>2) saturation_flag=1;
		for (int i=0; i<TDC_max-10; i++){
			if (adc_[i]>(0.4*(ADC_max_-ped_lv)+ped_lv)) {fardoublepeak_flag=1; break;}
		}
		for (int i=TDC_max+30; i<350; i++){
			if (adc_[i]>(0.4*(ADC_max_-ped_lv)+ped_lv)) {fardoublepeak_flag=1; break;}
		}
		for (int i=TDC_max-3; i<TDC_max; i++){
			if (adc_[i]<adc_[i-1]) {closedoublepeak_flag=1; break;}
		}
		for (int i=TDC_max+1; i<TDC_max+3; i++){
			if (adc_[i]<adc_[i+1]) {closedoublepeak_flag=1; break;}
		}
		
		if ((n_fit>85)&&(n_fit<100)&&(ped_lv<400.)&&(saturation_flag==0)&&(fardoublepeak_flag==0)&&(closedoublepeak_flag==0)){
		//if ((n_fit>85)&&(n_fit<100)&&(ped_lv<400.)&&(saturation_flag==0)&&(fardoublepeak_flag==0)){
			int n_non_fit=sampling_len_-n_fit;
			double ADC_fit[n_fit];
			double t_fit[n_fit];
			double dy[n_fit];//necessary to make sure chisquare is var/expt in TGraph fit, not var/1
			double ADC_non_fit[n_non_fit];
			double t_non_fit[n_non_fit];
			for (int i=0; i<TDC_start_peak; i++){
				ADC_non_fit[i]=(double)adc_[i]-ped_lv;
				t_non_fit[i]=PERTDC*(double)i;
			}
			for (int i=TDC_start_peak; i<(TDC_end_peak+1); i++){
				ADC_fit[i-TDC_start_peak]=(double)adc_[i]-ped_lv;
				t_fit[i-TDC_start_peak]=PERTDC*(double)i;
				dy[i-TDC_start_peak]=TMath::Sqrt(fabs(((double)adc_[i]-ped_lv))/singlePE)*singlePE;
				//dy[i-TDC_start_peak]=1.;
			}
			for (int i=(TDC_end_peak+1); i<sampling_len_; i++){
				ADC_non_fit[i-n_fit]=(double)adc_[i]-ped_lv;
				t_non_fit[i-n_fit]=PERTDC*(double)i;
			}
			TCanvas* c=new TCanvas(Form("Instance%i%s%i", index+1, labelname.c_str(), func_mode), "", 1600, 600);
			c->Divide(2,1);
			c->cd(1);
			gStyle->SetOptFit(1111);
			TGraphErrors* g1_fit=new TGraphErrors(n_fit, t_fit, ADC_fit, 0, dy);
			//TH1D* g1_fit=new TH1D(Form("fit_instance%i%s", index+1, labelname.c_str()), Form("Fitting for Instance #%i (%s); t [ns]; ADC", index+1, labelname.c_str()), n_fit, t_fit[0]-PERTDC/2., t_fit[n_fit-1]+PERTDC/2.);
			//for (int i=0; i<n_fit; i++){
				//g1_fit->Fill(t_fit[i], ADC_fit[i]);
			//}
			TGraph* g1=new TGraph(n_fit, t_fit, ADC_fit);
			TGraph* g2=new TGraph(n_non_fit, t_non_fit, ADC_non_fit);
			g1_fit->SetMarkerColor(kWhite);
			g1_fit->SetLineColor(kWhite);
			g1->SetMarkerStyle(24);
			g1->SetMarkerColor(kBlack);
			g1->SetLineColor(kBlack);
			g2->SetMarkerStyle(24);
			g2->SetMarkerColor(kBlue);
			g2->SetLineColor(kBlue);
			g1->SetMarkerSize(.7);
			g1->SetLineWidth(2);
			g1->SetFillColor(0);
			g2->SetMarkerSize(.7);
			g2->SetLineWidth(2);
			g2->SetFillColor(0);
			
			TMultiGraph* mg=new TMultiGraph();
			mg->SetTitle(Form("Fitting for Instance #%i (%s); t [ns]; ADC", index+1, labelname.c_str()));
			mg->Add(g1);
			mg->Add(g2);
			mg->Add(g1_fit);
			mg->Draw("AP");
			mg->GetYaxis()->SetTitleOffset(1.4);
			
			//fit g1 and extract fitting parameters, chi2, 50% rising edge
			//func_mode: 1=demg; 2=landau
			TF1* func1;
			if (func_mode==1){
				func1=new TF1(Form("Instance%i%sdemg", index+1, labelname.c_str()), "[0]/[2]*[1]*TMath::Exp(.5*[5]*[5]/[2]/[2]-(x-[4])/[2])*TMath::Erfc(1./TMath::Sqrt(2.)*([5]/[2]-(x-[4])/[5]))+[0]/[3]*(1.-[1])*TMath::Exp(.5*[5]*[5]/[3]/[3]-(x-[4])/[3])*TMath::Erfc(1./TMath::Sqrt(2.)*([5]/[3]-(x-[4])/[5]))", ((double)TDC_start_peak-0.5)*PERTDC, ((double)TDC_end_peak+0.5)*PERTDC);
				func1->SetParameter(0, 10.*(ADC_max_-ped_lv));
				func1->SetParameter(1, .5);
				func1->SetParameter(2, 15.);
				func1->SetParameter(3, 70.);
				func1->SetParameter(4, t_max_);
				func1->SetParameter(5, 2.5);
				
				func1->SetParLimits(0, 0.*(ADC_max_-ped_lv), 20.*(ADC_max_-ped_lv));
				func1->SetParLimits(1, 0., 1.);
				func1->SetParLimits(2, 1., 20.);
				func1->SetParLimits(3, 10., 150.);
				func1->SetParLimits(4, t_max_-10.*PERTDC, t_max_+10.*PERTDC);
				func1->SetParLimits(5, 0., 30.);
				
				func1->SetParName(0, "A");
				func1->SetParName(1, "r");
				func1->SetParName(2, "#tau_{1}");
				func1->SetParName(3, "#tau_{2}");
				func1->SetParName(4, "#mu");
				func1->SetParName(5, "#sigma");
			}
			else{
				func1=new TF1(Form("Instance%i%sland", index+1, labelname.c_str()), "[0]*TMath::Landau(x,[1],[2],1)", ((double)TDC_start_peak-0.5)*PERTDC, ((double)TDC_end_peak+0.5)*PERTDC);
				func1->SetParameter(0, ADC_max_-ped_lv);
				func1->SetParameter(1, t_max_);
				func1->SetParameter(2, 3.);
				
				func1->SetParLimits(0, 0.*(ADC_max_-ped_lv), 100000.*(ADC_max_-ped_lv));
				func1->SetParLimits(1, t_max_-20., t_max_+200.);
				func1->SetParLimits(2, 0.0001, 30.);
				
				func1->SetParName(0, "A");
				func1->SetParName(1, "#mu");
				func1->SetParName(2, "#sigma");
			}
			TFitResultPtr frp1=g1_fit->Fit(func1, "QRS");
			fp_A_0_=ped_lv;
			if (func_mode==1){
				fp_demg_A_=func1->GetParameter(0);
				fp_demg_r_=func1->GetParameter(1);
				fp_demg_tau1_=func1->GetParameter(2);
				fp_demg_tau2_=func1->GetParameter(3);
				fp_demg_mu_=func1->GetParameter(4);
				fp_demg_sigma_=func1->GetParameter(5);
				fp_demg_chi2_=func1->GetChisquare();
			}
			else {
				fp_land_A_=func1->GetParameter(0);
				fp_land_mu_=func1->GetParameter(1);
				fp_land_sigma_=func1->GetParameter(2);
				fp_land_chi2_=func1->GetChisquare();
			}
			t_50_=func1->GetX(func1->GetMaximum(((double)TDC_start_peak-0.5)*PERTDC, ((double)TDC_end_peak+0.5)*PERTDC)*0.5, t_max_-15.*PERTDC, t_max_);
			if (func_mode==1){
				if (((fp_demg_chi2_/ADC_max_)<CHI2THRESHOLD)&&(fp_demg_mu_>50.)&&(fabs(t_50_-t_50_ref_)<20.)){
					is_good_fit_demg_=1;
					if (plotflag==2) plotflag=0;
				}
				else {
					if (plotflag==2) plotflag=1;
				}
			}
			else {
				if (((fp_land_chi2_/ADC_max_)<CHI2THRESHOLD)&&(fp_land_mu_>50.)&&(fabs(t_50_-t_50_ref_)<20.)){
					is_good_fit_land_=1;
					if (plotflag==2) plotflag=0;
				}
				else {
					if (plotflag==2) plotflag=1;
				}
			}
			
			//if (plotflag) cout<<is_good_fit_<<" "<<endl;
			
			//plot if needed
			func1->SetLineColor(kRed);
			func1->DrawClone("SAME");
			TLegend* leg=new TLegend(.65, .50, .9, .35, "");
			leg->SetFillColor(0);
			leg->AddEntry(g2, "data points");
			leg->AddEntry(g1, "fitting range");
			leg->AddEntry(func1, "fit");
			leg->Draw("SAME");
			if (plotflag){
				c->cd(2);
				TPad *pad1 = new TPad("pad1","",0,0.4,1,1);
				pad1->SetFillColor(0);
				pad1->SetFillStyle(4000);
				pad1->SetFrameFillStyle(0);
				pad1->Draw();
				pad1->cd();
				TGraphErrors* g1_fit_dup=new TGraphErrors(n_fit, t_fit, ADC_fit, 0, dy);
				g1_fit_dup->SetMarkerColor(kBlack);
				g1_fit_dup->SetLineColor(kBlack);
				g1_fit_dup->SetTitle(";t [ns]; ADC");
				g1_fit_dup->Draw("AP");
				g1_fit_dup->GetXaxis()->SetLimits(PERTDC*((double)TDC_start_peak-0.5),PERTDC*((double)TDC_end_peak+0.5));
				g1_fit_dup->GetXaxis()->SetLabelSize(0.075);
				g1_fit_dup->GetYaxis()->SetLabelSize(0.075);
				g1_fit_dup->GetXaxis()->SetTitleSize(0.075);
				g1_fit_dup->GetYaxis()->SetTitleSize(0.075);
				g1_fit_dup->GetXaxis()->SetTitleOffset(0.6);
				g1_fit_dup->GetYaxis()->SetTitleOffset(0.6);
				func1->DrawClone("SAME");
				c->cd(2);
				TPad *pad2 = new TPad("pad2","",0,0,1,0.4);
				pad2->SetFillColor(0);
				pad2->SetFillStyle(4000);
				pad2->SetFrameFillStyle(0);
				pad2->Draw();
				pad2->cd();
				double error[n_fit];
				double xline[n_fit];
				for (int i=0; i<n_fit; i++){
					error[i]=(ADC_fit[i]-func1->Eval(t_fit[i]))/dy[i];
					xline[i]=PERTDC/2.;
				}
				TGraphErrors* err=new TGraphErrors(n_fit, t_fit, error, xline, 0);
				err->SetTitle(";;[ADC(t)-f(t)]/#sigma");
				err->SetMarkerColor(kBlue);
				err->SetLineColor(kBlue);
				err->Draw("AP");
				err->GetXaxis()->SetLimits(PERTDC*((double)TDC_start_peak-0.5),PERTDC*((double)TDC_end_peak+0.5));
				err->GetXaxis()->SetLabelSize(0.125);
				err->GetYaxis()->SetLabelSize(0.125);
				err->GetXaxis()->SetTitleSize(0.125);
				err->GetYaxis()->SetTitleSize(0.125);
				err->GetYaxis()->SetTitleOffset(0.4);
				TLine* line=new TLine(PERTDC*((double)TDC_start_peak-0.5),0.,PERTDC*((double)TDC_end_peak+0.5),0.);
				line->SetLineWidth(2);
				line->SetLineColor(kBlack);
				line->Draw("SAME");
				c->Print(pdfout.c_str(), "pdf");
			}
			else{
				delete leg;
				//delete func1;
				delete g1;
				delete g2;
				delete mg;
				//delete g1_fit;
				delete c;
			}
			
			return;
		}
		else{
			if (plotflag!=0){
				double ADC_all[sampling_len_];
				double time_all[sampling_len_];
				for (int i=0; i<sampling_len_; i++){
					ADC_all[i]=(double)adc_[i];
					time_all[i]=PERTDC*(double)i;
				}
				TCanvas* c=new TCanvas(Form("Instance%i%s%i", index+1, labelname.c_str(),func_mode), "", 1600, 600);
				c->Divide(2,1);
				c->cd(1);
				TGraph* g1=new TGraph(sampling_len_, time_all, ADC_all);
				g1->SetTitle(Form("Fitting for Instance #%i (%s) range_TDC(%i,%i) ped %f n_fit=%i %i %i %i; t [ns]; ADC", index+1, labelname.c_str(),TDC_start_peak,TDC_end_peak,ped_lv, n_fit, saturation_flag, fardoublepeak_flag, closedoublepeak_flag));
				g1->SetMarkerStyle(24);
				g1->SetMarkerColor(kBlack);
				g1->SetLineColor(kBlack);
				g1->SetMarkerSize(.7);
				g1->SetLineWidth(2);
				g1->SetFillColor(0);
				g1->Draw("ALP");
				g1->GetYaxis()->SetTitleOffset(1.4);
				c->Print(pdfout.c_str(), "pdf");
			}
			return;
		}
	}
};

class event{
private:
	channel near_;
	channel far_;
	int near_is_good_;
	int far_is_good_;
	double t_diff_max_;
	double t_diff_50_;
public:
	void init(){
		near_.init();
		far_.init();
		near_is_good_=0;
		far_is_good_=0;
		t_diff_max_=0.;
		t_diff_50_=0.;
		return;
	}
	event(){
		near_=channel();
		far_=channel();
		init();
	}
	channel near() const {return near_;}
	channel far() const {return far_;}
	int near_is_good(){return near_is_good_;}
	int far_is_good(){return far_is_good_;}
	double t_diff_max(){return t_diff_max_;}
	double t_diff_50(){return t_diff_50_;}
	
	void setProcessedBranches(TTree* processed){/////coarse
		near_.setProcessedBranches(processed, 'n');
		far_.setProcessedBranches(processed, 'f');
		processed->Branch("t_diff_max", &t_diff_max_, "t_diff_max/D");
		processed->Branch("t_diff_50", &t_diff_50_, "t_diff_50/D");
		return;
	}
	void setRawBranchAddresses(TTree* raw){
		near_.setRawBranchAddresses(raw, 'n');
		far_.setRawBranchAddresses(raw, 'f');
		raw->SetBranchAddress("n_quality", &near_is_good_);
		raw->SetBranchAddress("f_quality", &far_is_good_);
		return;
	}
	void setCoarseBranchAddresses(TTree* coarse){
		near_.setCoarseBranchAddresses(coarse, 'n');
		far_.setCoarseBranchAddresses(coarse, 'f');
		return;
	}
	void setProcessedBranchAddresses(TTree* processed){
		near_.setProcessedBranchAddresses(processed, 'n');
		far_.setProcessedBranchAddresses(processed, 'f');
		processed->SetBranchAddress("t_diff_max", &t_diff_max_);
		processed->SetBranchAddress("t_diff_50", &t_diff_50_);
		return;
	}
	int process(int index, int func_mode, string const &pdfout, int plotflag=0){//return: 0,both good fit; 1,bad fit; 2,bad data
		if ((near_is_good_)&&(far_is_good_)){
			near_.fit_func(func_mode, index, 'n', pdfout, plotflag);
			far_.fit_func(func_mode, index, 'f', pdfout, plotflag);
			int near_good_fit, far_good_fit;
			if (func_mode==1){
				near_good_fit=near_.is_good_fit_demg();
				far_good_fit=far_.is_good_fit_demg();
			}
			else {
				near_good_fit=near_.is_good_fit_land();
				far_good_fit=far_.is_good_fit_land();
			}
			if (near_good_fit&&far_good_fit){
				t_diff_max_=far_.t_max()-near_.t_max();
				t_diff_50_=far_.t_50()-near_.t_50();
				return 0;
			}
			else return 1;
		}
		else return 2;
	}
	void show_both(int index, string const &pdfout){
		int sampling_length=N_width*8;
		double time[sampling_length];
		double adc_n[sampling_length];
		double adc_f[sampling_length];
		for (int i=0; i<sampling_length; i++){
			time[i]=(double)i*PERTDC;
			adc_n[i]=(double)near_.adc_i(i);
			adc_f[i]=(double)far_.adc_i(i);
		}
		bichannel_plot(pdfout, adc_n, adc_f, time, sampling_length, index);
		return;
	}
};
