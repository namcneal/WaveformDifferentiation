# include "class.cpp"

using namespace std;

int main(int argc, char **argv){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;

	string name (argv[1]);
	//string plot_name=name+"_plot_sample_fit.pdf";
	string plot_name1=name+"_plot_sample_bad_fit_demg.pdf";
	string plot_name2=name+"_plot_sample_bad_fit_land.pdf";
	TCanvas* c0=new TCanvas("0", "", 1600, 600);
	c0->Print((plot_name1+"[").c_str(),"pdf");

	TFile* new_file1=new TFile((name+"_demg_processed_data.root").c_str(), "recreate");

	new_file1->cd();
	TTree* new_tree1=new TTree((name+"_demg_processed").c_str(),(name+"_demg_processed").c_str());
	TFile* raw_file=new TFile((name+"_raw_data.root").c_str(), "READ");
	TTree* raw_tree=(TTree*)raw_file->Get(Form("%s_raw", name.c_str()));
	TFile* coarse_file=new TFile((name+"_coarse_time.root").c_str(), "READ");
	TTree* coarse_tree=(TTree*)coarse_file->Get(Form("%s_coarse_time", name.c_str()));
	event* dynamic1=new event();
	dynamic1->setProcessedBranches(new_tree1);
	dynamic1->setRawBranchAddresses(raw_tree);
	dynamic1->setCoarseBranchAddresses(coarse_tree);
	int n_tot=raw_tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<name<<"_raw_data.root"<<endl;

	cout<<"Start fitting with Double Exponentially Modified Gaussian..."<<endl;

	//int flag=0;
	int n_bad_data=0;
	int n_bad_fit=0;
	int n_fit=0;

	TH1D* h_chi2n1=new TH1D("h_chi2n1", "Distribution of #chi^{2} in Asymmetric Gaussian Fit; #chi^{2}; count", 100, 0., 180.*20.*20.);
	TH1D* h_chi2f1=new TH1D("h_chi2f1", "Distribution of #chi^{2} in Asymmetric Gaussian Fit; #chi^{2}; count", 100, 0., 180.*20.*20.);

	TH1D* h_t_diff_501=new TH1D("h_t_diff_501", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);
	TH1D* h_t_diff_max1=new TH1D("h_t_diff_max1", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);

	//for (int i=0; i<n_tot; i++){
	for (int i=0; i<500; i++){
		int plotflag=0;
		if (i<50) plotflag=1;
		else plotflag=2;
		dynamic1->init();
		raw_tree->GetEntry(i);
		coarse_tree->GetEntry(i);
		int flag=dynamic1->process(i, 1, plot_name1, plotflag);
		new_tree1->Fill();
		if (flag==2) n_bad_data++;
		else{
			//h_chi2n1->Fill((dynamic1->near()).fp_demg_chi2());
			//h_chi2f1->Fill((dynamic1->far()).fp_demg_chi2());
			if (flag==0){
				n_fit++;
				h_chi2n1->Fill((dynamic1->near()).fp_demg_chi2());
				h_chi2f1->Fill((dynamic1->far()).fp_demg_chi2());
				h_t_diff_501->Fill(dynamic1->t_diff_50());
				h_t_diff_max1->Fill(dynamic1->t_diff_max());
			}
			else n_bad_fit++;
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_bad_data<<" events have empty channel."<<endl;
	cout<<"    "<<n_bad_fit<<" events have bad fitting in either channel."<<endl;
	cout<<"    "<<n_fit<<" events have proper fits and time difference of signals were calculated."<<endl;
	c0->Print((plot_name1+"]").c_str(),"pdf");

	new_file1->cd();
	new_tree1->Write();
	delete dynamic1;

	TFile* new_file2=new TFile((name+"_land_processed_data.root").c_str(), "recreate");
	TTree* new_tree2=new TTree((name+"_land_processed").c_str(),(name+"_land_processed").c_str());
	new_file2->cd();
	event* dynamic2=new event();
	dynamic2->setProcessedBranches(new_tree2);
	dynamic2->setRawBranchAddresses(raw_tree);
	dynamic2->setCoarseBranchAddresses(coarse_tree);
	cout<<"Start fitting with Landau Function..."<<endl;

	//int flag=0;
	n_bad_data=0;
	n_bad_fit=0;
	n_fit=0;
	c0->Print((plot_name2+"[").c_str(),"pdf");

	TH1D* h_chi2n2=new TH1D("h_chi2n2", "Distribution of #chi^{2} in Asymmetric Gaussian Fit; #chi^{2}; count", 100, 0., 180.*20.*20.);
	TH1D* h_chi2f2=new TH1D("h_chi2f2", "Distribution of #chi^{2} in Asymmetric Gaussian Fit; #chi^{2}; count", 100, 0., 180.*20.*20.);

	TH1D* h_t_diff_502=new TH1D("h_t_diff_502", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);
	TH1D* h_t_diff_max2=new TH1D("h_t_diff_max2", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);

	//for (int i=0; i<n_tot; i++){
	for (int i=0; i<500; i++){
		int plotflag=0;
		if (i<50) plotflag=1;
		else plotflag=2;
		dynamic2->init();
		raw_tree->GetEntry(i);
		coarse_tree->GetEntry(i);
		int flag=dynamic2->process(i, 2, plot_name2, plotflag);
		new_tree2->Fill();
		if (flag==2) n_bad_data++;
		else{
			//h_chi2n2->Fill((dynamic2->near()).fp_land_chi2());
			//h_chi2f2->Fill((dynamic2->far()).fp_land_chi2());
			if (flag==0){
				n_fit++;
				h_chi2n2->Fill((dynamic2->near()).fp_land_chi2());
				h_chi2f2->Fill((dynamic2->far()).fp_land_chi2());
				h_t_diff_502->Fill(dynamic2->t_diff_50());
				h_t_diff_max2->Fill(dynamic2->t_diff_max());
			}
			else n_bad_fit++;
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_bad_data<<" events have empty channel."<<endl;
	cout<<"    "<<n_bad_fit<<" events have bad fitting in either channel."<<endl;
	cout<<"    "<<n_fit<<" events have proper fits and time difference of signals were calculated."<<endl;
	c0->Print((plot_name2+"]").c_str(),"pdf");

	new_file2->cd();
	new_tree2->Write();
	delete dynamic2;

	TCanvas* c11=new TCanvas("11", "", 800, 600);
	c11->cd();
	c11->SetLogy();
	gStyle->SetOptStat(110010);
	h_chi2n1->SetLineColor(kRed);
	h_chi2f1->SetLineColor(kBlue);
	h_chi2n1->Draw();
	c11->Update();
	TLegend* leg111=new TLegend(.6, .6, .9, .4, "");
	leg111->SetFillColor(0);
	leg111->AddEntry(h_chi2n1, "near channel");
	leg111->Draw("SAME");
	c11->Print((name+"_demg_chi2_dist.pdf(").c_str(),"pdf");
	h_chi2f1->Draw();
	c11->Update();
	TLegend* leg121=new TLegend(.6, .6, .9, .4, "");
	leg121->SetFillColor(0);
	leg121->AddEntry(h_chi2f1, "far  channel");
	leg121->Draw("SAME");
	c11->Print((name+"_demg_chi2_dist.pdf").c_str(),"pdf");
	gStyle->SetOptStat(0);
	h_chi2n1->Draw();
	c11->Update();
	h_chi2f1->Draw("SAME");
	c11->Update();
	TLegend* leg131=new TLegend(.6, .9, .9, .7, "");
	leg131->SetFillColor(0);
	leg131->AddEntry(h_chi2n1, "near channel");
	leg131->AddEntry(h_chi2f1, "far  channel");
	leg131->Draw("SAME");
	c11->Print((name+"_demg_chi2_dist.pdf)").c_str(),"pdf");

	TCanvas* c21=new TCanvas("21", "", 800, 600);
	c21->cd();
	gStyle->SetOptStat(110010);
	h_t_diff_501->SetLineColor(kRed);
	h_t_diff_max1->SetLineColor(kBlack);
	h_t_diff_max1->SetMinimum(0.8);
	h_t_diff_max1->Draw();
	h_t_diff_501->Draw("SAME");
	TLegend* leg21=new TLegend(.6, .9, .9, .7, "");
	leg21->SetFillColor(0);
	leg21->AddEntry(h_t_diff_501, "50% rising edge");
	leg21->AddEntry(h_t_diff_max1, "maximum point (discrete)");
	leg21->Draw("SAME");
	c21->Print((name+"_demg_time_dist.pdf(").c_str(),"pdf");
	c21->SetLogy();
	c21->Update();
	c21->Print((name+"_demg_time_dist.pdf)").c_str(),"pdf");

	TCanvas* c12=new TCanvas("12", "", 800, 600);
	c12->cd();
	c12->SetLogy();
	gStyle->SetOptStat(110010);
	h_chi2n2->SetLineColor(kRed);
	h_chi2f2->SetLineColor(kBlue);
	h_chi2n2->Draw();
	c12->Update();
	TLegend* leg112=new TLegend(.6, .6, .9, .4, "");
	leg112->SetFillColor(0);
	leg112->AddEntry(h_chi2n2, "near channel");
	leg112->Draw("SAME");
	c12->Print((name+"_land_chi2_dist.pdf(").c_str(),"pdf");
	h_chi2f2->Draw();
	c12->Update();
	TLegend* leg122=new TLegend(.6, .6, .9, .4, "");
	leg122->SetFillColor(0);
	leg122->AddEntry(h_chi2f2, "far  channel");
	leg122->Draw("SAME");
	c12->Print((name+"_land_chi2_dist.pdf").c_str(),"pdf");
	gStyle->SetOptStat(0);
	h_chi2n2->Draw();
	c12->Update();
	h_chi2f2->Draw("SAME");
	c12->Update();
	TLegend* leg132=new TLegend(.6, .9, .9, .7, "");
	leg132->SetFillColor(0);
	leg132->AddEntry(h_chi2n2, "near channel");
	leg132->AddEntry(h_chi2f2, "far  channel");
	leg132->Draw("SAME");
	c12->Print((name+"_land_chi2_dist.pdf)").c_str(),"pdf");

	TCanvas* c22=new TCanvas("22", "", 800, 600);
	c22->cd();
	gStyle->SetOptStat(110010);
	h_t_diff_502->SetLineColor(kRed);
	h_t_diff_max2->SetLineColor(kBlack);
	h_t_diff_max2->SetMinimum(0.8);
	h_t_diff_max2->Draw();
	h_t_diff_502->Draw("SAME");
	TLegend* leg22=new TLegend(.6, .9, .9, .7, "");
	leg22->SetFillColor(0);
	leg22->AddEntry(h_t_diff_502, "50% rising edge");
	leg22->AddEntry(h_t_diff_max2, "maximum point (discrete)");
	leg22->Draw("SAME");
	c22->Print((name+"_land_time_dist.pdf(").c_str(),"pdf");
	c22->SetLogy();
	c22->Update();
	c22->Print((name+"_land_time_dist.pdf)").c_str(),"pdf");

	raw_file->Close();
	new_file1->Close();
	new_file2->Close();

	return 0;
}
