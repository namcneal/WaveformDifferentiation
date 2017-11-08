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
	TFile* processed=new TFile((name+"_land_processed_data.root").c_str(), "READ");
	TTree* tree=(TTree*)processed->Get(Form("%s_land_processed", name.c_str()));

	string plot_name=name+"_land_parameters.pdf";
	TCanvas* c0=new TCanvas("0", "", 800, 600);
	c0->Print((plot_name+"[").c_str(),"pdf");
	
	event* dynamic=new event();
	dynamic->setProcessedBranchAddresses(tree);
	int n_tot=tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<name<<"_land_processed_data.root"<<endl;

	//gStyle->SetOptStat(111111);
	
	TH1D* h_t_diff_50=new TH1D("h_t_diff_50", "Distribution of Time Difference; #Delta t  [ns]; count", 140, -80., 200.);
	TH1D* h_t_diff_max=new TH1D("h_t_diff_max", "Distribution of Time Difference; #Delta t  [ns]; count", 140, -80., 200.);
	TH1D* h_fp_sigma=new TH1D("h_fp_sigma","Distribution of #sigma in fit; A [adc]; count", 100, 2., 7.);
	TH1D* h_fp_sigma_n=new TH1D("h_fp_sigma_n","Distribution of #sigma in fit; A [adc]; count", 100, 2., 7.);
	TH1D* h_fp_sigma_p=new TH1D("h_fp_sigma_p","Distribution of #sigma in fit; A [adc]; count", 100, 2., 7.);
	
	TH2D* h_t_diff_50_sigma=new TH2D("h_t_diff_50_sigma", "Correlation of #sigma vs. #Delta t (50% rising edge); #Delta t [ns]; #sigma [ns]", 140, -80., 200., 100, 2., 7.);
	
	TH2D* h_ADC_max_fp_chi2=new TH2D("h_ADC_max_fp_chi2", "Correlation of #chi^{2} vs. ADC_max ; ADC_max; #chi^{2}", 100, 0., 10000., 100, 0., 200*20.*20.);

	for (int i=0; i<n_tot; i++){
		dynamic->init();
		tree->GetEntry(i);
		if (((dynamic->near()).is_good_fit_land())&&((dynamic->far()).is_good_fit_land())){
			double tdiff50=dynamic->t_diff_50();
			h_t_diff_50->Fill(tdiff50);
			h_t_diff_max->Fill(dynamic->t_diff_max());
			double fp_sigma_far=(dynamic->far()).fp_land_sigma();
			h_t_diff_50_sigma->Fill(tdiff50, fp_sigma_far);
			h_fp_sigma->Fill(fp_sigma_far);
			h_ADC_max_fp_chi2->Fill((dynamic->near()).ADC_max(),(dynamic->near()).fp_land_chi2());
			h_ADC_max_fp_chi2->Fill((dynamic->far()).ADC_max(),(dynamic->far()).fp_land_chi2());
			if ((tdiff50<15)&&(tdiff50>-20)){
				h_fp_sigma_p->Fill(fp_sigma_far);
			}
			if (tdiff50>30){
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
	h_fp_sigma_n->SetLineColor(kRed);
	h_fp_sigma_p->SetLineColor(kBlue);
	h_fp_sigma->SetLineColor(kBlack);
	h_fp_sigma->Draw();
	h_fp_sigma_n->Draw("SAME");
	h_fp_sigma_p->Draw("SAME");
	h_fp_sigma->Draw("SAMES");
	TLegend* leg1=new TLegend(.75, .9, .95, .75, "");
	leg1->SetFillColor(0);
	leg1->AddEntry(h_fp_sigma_n, "neutron");
	leg1->AddEntry(h_fp_sigma_p, "photon");
	leg1->AddEntry(h_fp_sigma, "sum");
	leg1->Draw("SAME");
	c1->Print(plot_name.c_str(),"pdf");
	c1->SetLogy();
	c1->Print(plot_name.c_str(),"pdf");
	
	c1->SetLogy(0);
	TCanvas* c2=new TCanvas("2", "", 800, 600);
	c2->cd();
	gStyle->SetOptStat(0);
	h_t_diff_50_sigma->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	h_ADC_max_fp_chi2->Draw("COLZ");
	c2->Print(plot_name.c_str(),"pdf");
	
	c0->Print((plot_name+"]").c_str(),"pdf");

	processed->Close();
	return 0;
}
