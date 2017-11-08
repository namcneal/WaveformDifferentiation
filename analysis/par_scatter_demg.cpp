# include "class.cpp"

using namespace std;

///////////////////////////////////////////

int main(int argc, char **argv){
	gROOT->Reset();
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
	
	string name (argv[1]);
	TFile* processed=new TFile((name+"_demg_processed_data.root").c_str(), "READ");
	TTree* tree=(TTree*)processed->Get(Form("%s_demg_processed", name.c_str()));

	string plot_name=name+"_demg_parameters.pdf";
	TCanvas* c0=new TCanvas("0", "", 800, 600);
	c0->Print((plot_name+"[").c_str(),"pdf");
	
	event* dynamic=new event();
	dynamic->setProcessedBranchAddresses(tree);
	int n_tot=tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<name<<"_demg_processed_data.root"<<endl;

	//gStyle->SetOptStat(111111);
	
	TH1D* h_t_diff_50=new TH1D("h_t_diff_50", "Distribution of Time Difference; #Delta t  [ns]; count", 140, -80., 200.);
	TH1D* h_t_diff_max=new TH1D("h_t_diff_max", "Distribution of Time Difference; #Delta t  [ns]; count", 140, -80., 200.);
	
	TH1D* h_fp_A=new TH1D("h_fp_A","Distribution of amplitude in fit; A [adc]; count", 100, 0., 100000.);
	TH1D* h_fp_r=new TH1D("h_fp_r","Distribution of r in fit; r; count", 100, 0., 1.);
	TH1D* h_fp_tau1=new TH1D("h_fp_tau1","Distribution of #tau_{1} in fit; #tau_{1} [ns]; count", 100, 1., 20.);
	TH1D* h_fp_tau2=new TH1D("h_fp_tau2","Distribution of #tau_{2} in fit; #tau_{2} [ns]; count", 100, 10., 150.);
	TH1D* h_fp_sigma=new TH1D("h_fp_sigma","Distribution of #sigma in fit; #sigma [ns]; count", 100, 0., 5.);
	
	TH1D* h_fp_A_n=new TH1D("h_fp_A_n","Distribution of amplitude in fit; A [adc]; count", 100, 0., 100000.);
	TH1D* h_fp_r_n=new TH1D("h_fp_r_n","Distribution of r in fit; r; count", 100, 0., 1.);
	TH1D* h_fp_tau1_n=new TH1D("h_fp_tau1_n","Distribution of #tau_{1} in fit; #tau_{1} [ns]; count", 100, 1., 20.);
	TH1D* h_fp_tau2_n=new TH1D("h_fp_tau2_n","Distribution of #tau_{2} in fit; #tau_{2} [ns]; count", 100, 10., 150.);
	TH1D* h_fp_sigma_n=new TH1D("h_fp_sigma_n","Distribution of #sigma in fit; #sigma [ns]; count", 100, 0., 5.);
	
	TH1D* h_fp_A_p=new TH1D("h_fp_A_p","Distribution of amplitude in fit; A [adc]; count", 100, 0., 100000.);
	TH1D* h_fp_r_p=new TH1D("h_fp_r_p","Distribution of r in fit; r; count", 100, 0., 1.);
	TH1D* h_fp_tau1_p=new TH1D("h_fp_tau1_p","Distribution of #tau_{1} in fit; #tau_{1} [ns]; count", 100, 1., 20.);
	TH1D* h_fp_tau2_p=new TH1D("h_fp_tau2_p","Distribution of #tau_{2} in fit; #tau_{2} [ns]; count", 100, 10., 150.);
	TH1D* h_fp_sigma_p=new TH1D("h_fp_sigma_p","Distribution of #sigma in fit; #sigma [ns]; count", 100, 0., 5.);
	
	TH2D* h_t_diff_50_A=new TH2D("h_t_diff_50_A", "Correlation of A vs. #Delta t (50% rising edge); #Delta t [ns]; A", 140, -80., 200., 100, 0., 100000.);
	TH2D* h_t_diff_50_r=new TH2D("h_t_diff_50_r", "Correlation of r vs. #Delta t (50% rising edge); #Delta t [ns]; r", 140, -80., 200., 100, 0., 1.);
	TH2D* h_t_diff_50_tau1=new TH2D("h_t_diff_50_tau1", "Correlation of #tau_{1} vs. #Delta t (50% rising edge); #Delta t [ns]; #tau_{1} [ns]", 140, -80., 200., 100, 1., 20.);
	TH2D* h_t_diff_50_tau2=new TH2D("h_t_diff_50_tau2", "Correlation of #tau_{2} vs. #Delta t (50% rising edge); #Delta t [ns]; #tau_{2} [ns]", 140, -80., 200., 100, 10., 150.);
	TH2D* h_t_diff_50_sigma=new TH2D("h_t_diff_50_sigma", "Correlation of #sigma vs. #Delta t (50% rising edge); #Delta t [ns]; #sigma [ns]", 140, -80., 200., 100, 0., 5.);

	TH2D* h_ADC_max_fp_chi2=new TH2D("h_ADC_max_fp_chi2", "Correlation of #chi^{2} vs. ADC_max ; ADC_max; #chi^{2}", 100, 0., 100000., 100, 0., 200*20.*20.);

	for (int i=0; i<n_tot; i++){
		dynamic->init();
		tree->GetEntry(i);
		if (((dynamic->near()).is_good_fit_demg())&&((dynamic->far()).is_good_fit_demg())){
			double tdiff50=dynamic->t_diff_50();
			h_t_diff_50->Fill(tdiff50);
			h_t_diff_max->Fill(dynamic->t_diff_max());
			double fp_A_far=(dynamic->far()).fp_demg_A();
			double fp_r_far=(dynamic->far()).fp_demg_r();
			double fp_tau1_far=(dynamic->far()).fp_demg_tau1();
			double fp_tau2_far=(dynamic->far()).fp_demg_tau2();
			double fp_sigma_far=(dynamic->far()).fp_demg_sigma();
			h_fp_A->Fill(fp_A_far);
			h_fp_r->Fill(fp_r_far);
			h_fp_tau1->Fill(fp_tau1_far);
			h_fp_tau2->Fill(fp_tau2_far);
			h_fp_sigma->Fill(fp_sigma_far);
			h_t_diff_50_A->Fill(tdiff50, fp_A_far);
			h_t_diff_50_r->Fill(tdiff50, fp_r_far);
			h_t_diff_50_tau1->Fill(tdiff50, fp_tau1_far);
			h_t_diff_50_tau2->Fill(tdiff50, fp_tau2_far);
			h_t_diff_50_sigma->Fill(tdiff50, fp_sigma_far);
			h_ADC_max_fp_chi2->Fill((dynamic->near()).ADC_max(),(dynamic->near()).fp_demg_chi2());
			h_ADC_max_fp_chi2->Fill((dynamic->far()).ADC_max(),(dynamic->far()).fp_demg_chi2());
			if ((tdiff50<15)&&(tdiff50>-20)){
				h_fp_A_p->Fill(fp_A_far);
				h_fp_r_p->Fill(fp_r_far);
				h_fp_tau1_p->Fill(fp_tau1_far);
				h_fp_tau2_p->Fill(fp_tau2_far);
				h_fp_sigma_p->Fill(fp_sigma_far);
			}
			if (tdiff50>30){
				h_fp_A_n->Fill(fp_A_far);
				h_fp_r_n->Fill(fp_r_far);
				h_fp_tau1_n->Fill(fp_tau1_far);
				h_fp_tau2_n->Fill(fp_tau2_far);
				h_fp_sigma_n->Fill(fp_sigma_far);
			}
		}
	}

	TCanvas* c1=new TCanvas("1", "", 800, 600);
	c1->cd();
	//gStyle->SetOptStat(110011);
	h_t_diff_50->SetLineColor(kRed);
	h_t_diff_max->SetLineColor(kBlack);
	h_t_diff_max->SetMinimum(0.8);
	h_t_diff_max->Draw();
	h_t_diff_50->Draw("SAME");
	TLegend* leg0=new TLegend(.5, .88, .75, .74, "");
	leg0->SetFillColor(0);
	leg0->AddEntry(h_t_diff_50, "50% rising edge");
	leg0->AddEntry(h_t_diff_max, "maximum point (discrete)");
	leg0->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	c1->SetLogy(0);
	//gStyle->SetOptStat(111111);
	h_fp_A_n->SetLineColor(kRed);
	h_fp_A_p->SetLineColor(kBlue);
	h_fp_A->SetLineColor(kBlack);
	h_fp_A->Draw();
	h_fp_A_n->Draw("SAME");
	h_fp_A_p->Draw("SAME");
	h_fp_A->Draw("SAMES");
	TLegend* leg1=new TLegend(.75, .9, .95, .75, "");
	leg1->SetFillColor(0);
	leg1->AddEntry(h_fp_A_n, "neutron");
	leg1->AddEntry(h_fp_A_p, "photon");
	leg1->AddEntry(h_fp_A, "sum");
	leg1->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	c1->SetLogy(0);
	//gStyle->SetOptStat(111111);
	h_fp_r_n->SetLineColor(kRed);
	h_fp_r_p->SetLineColor(kBlue);
	h_fp_r->SetLineColor(kBlack);
	h_fp_r->Draw();
	h_fp_r_n->Draw("SAME");
	h_fp_r_p->Draw("SAME");
	h_fp_r->Draw("SAMES");
	TLegend* leg2=new TLegend(.75, .9, .95, .75, "");
	leg2->SetFillColor(0);
	leg2->AddEntry(h_fp_r_n, "neutron");
	leg2->AddEntry(h_fp_r_p, "photon");
	leg2->AddEntry(h_fp_r, "sum");
	leg2->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	c1->SetLogy(0);
	//gStyle->SetOptStat(111111);
	h_fp_tau1_n->SetLineColor(kRed);
	h_fp_tau1_p->SetLineColor(kBlue);
	h_fp_tau1->SetLineColor(kBlack);
	h_fp_tau1->Draw();
	h_fp_tau1_n->Draw("SAME");
	h_fp_tau1_p->Draw("SAME");
	h_fp_tau1->Draw("SAMES");
	TLegend* leg3=new TLegend(.75, .9, .95, .75, "");
	leg3->SetFillColor(0);
	leg3->AddEntry(h_fp_tau1_n, "neutron");
	leg3->AddEntry(h_fp_tau1_p, "photon");
	leg3->AddEntry(h_fp_tau1, "sum");
	leg3->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	c1->SetLogy(0);
	//gStyle->SetOptStat(111111);
	h_fp_tau2_n->SetLineColor(kRed);
	h_fp_tau2_p->SetLineColor(kBlue);
	h_fp_tau2->SetLineColor(kBlack);
	h_fp_tau2->Draw();
	h_fp_tau2_n->Draw("SAME");
	h_fp_tau2_p->Draw("SAME");
	h_fp_tau2->Draw("SAMES");
	TLegend* leg4=new TLegend(.75, .9, .95, .75, "");
	leg4->SetFillColor(0);
	leg4->AddEntry(h_fp_tau2_n, "neutron");
	leg4->AddEntry(h_fp_tau2_p, "photon");
	leg4->AddEntry(h_fp_tau2, "sum");
	leg4->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	c1->SetLogy(0);
	//gStyle->SetOptStat(111111);
	h_fp_sigma_n->SetLineColor(kRed);
	h_fp_sigma_p->SetLineColor(kBlue);
	h_fp_sigma->SetLineColor(kBlack);
	h_fp_sigma->Draw();
	h_fp_sigma_n->Draw("SAME");
	h_fp_sigma_p->Draw("SAME");
	h_fp_sigma->Draw("SAMES");
	TLegend* leg5=new TLegend(.75, .9, .95, .75, "");
	leg5->SetFillColor(0);
	leg5->AddEntry(h_fp_sigma_n, "neutron");
	leg5->AddEntry(h_fp_sigma_p, "photon");
	leg5->AddEntry(h_fp_sigma, "sum");
	leg5->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	TCanvas* c2=new TCanvas("2", "", 800, 600);
	c2->cd();
	gStyle->SetOptStat(0);
	h_t_diff_50_A->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	h_t_diff_50_r->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	h_t_diff_50_tau1->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	h_t_diff_50_tau2->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	h_t_diff_50_sigma->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	h_ADC_max_fp_chi2->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	
	c0->Print((plot_name+"]").c_str(),"pdf");

	processed->Close();
	return 0;
}
