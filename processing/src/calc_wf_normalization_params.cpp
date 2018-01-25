
#include "../include/calc_wf_normalization_params.h"
using namespace std;

/* Takes a raw ROOT data file and calculates the parameters used to normalize the waveforms. 
 * Stores the parameters (pedestals, amplitudes, max time, 50% rising edge) for each event in a new 
 * ROOT file called "normalization_params"
 */

// Read in the path for the ROOT files and the run name.
// Example of two arguments: ../../data take103
int calc_wf_normalization_params(string filePathPrefix, string fileName){
	// Define the name of the ROOT file with the raw data in it
    string raw_file_name = filePathPrefix + '/' + fileName + "_raw_data.root",
    	   raw_data_tree_name = fileName + "_raw";
	
	// Open the ROOT file and the Tree it contains with the raw data 
	TFile* raw_file=new TFile(raw_file_name.c_str(), "READ");
	TTree* raw_tree = (TTree*)raw_file->Get(raw_data_tree_name.c_str());
	
	// Set the name of the file and ROOT Tree that will contain the normalization parameters 
	string param_file_name = filePathPrefix + '/' + fileName + "_normalization_params.root",
		   new_data_tree_name = fileName + "_normalization_params";
		   
	// Actually create the ROOT file and ROOT Tree to contain the normalization parameters
	TFile* param_file = new TFile(param_file_name.c_str(), "recreate");
	TTree* param_tree = new TTree(new_data_tree_name.c_str(),(fileName+"_normalization_params").c_str());
	
	// Check to see if the raw data tree and normalization parameter tree were correctly opened
	if (!raw_tree){
		cout << "\n\nCannot find raw data tree" << endl << endl;
	}
	if (!param_tree){
		cout << "Could not creat the new " << fileName << "_course_time tree." << endl;
	}
	
	// Define variables into which we will read raw data
	bool is_good_near = 0, 
		is_good_far = 0;
		
	int adc_near[N_width*32] = {0},
		adc_far[N_width*32] = {0};
		
	// Set the connection between the ROOT raw data tree and the variables defined above 
	// through the variables' addresses in memory
	cout << "Setting branch addresses for the raw tree" << "...adc n";
	raw_tree->SetBranchAddress("adc_near", adc_near);
	cout << "...adc f";
	raw_tree->SetBranchAddress("adc_far", adc_far);
	cout << "...n quality";
	raw_tree->SetBranchAddress("near_quality", &is_good_near);
	cout << "...f quality." << endl;
	raw_tree->SetBranchAddress("far_quality", &is_good_far);
	
	// Check the total number of entries, storing the number 
	cout << "Geting the total number of entries." << endl;
	int n_tot = raw_tree->GetEntries();
	cout << "Read "<< n_tot <<" events from " << fileName << "_raw_data.root" << endl;
	
	// Variables to hold the data that will be moved into the normalization parameter file
	double ped_near = 0, 
	   ped_far = 0, 
	   amp_near = 0, 
	   amp_far = 0, 
	   t_max_near = 0, 
	   t_max_far = 0, 
	   t_50_near = 0, 
	   t_50_far = 0;

	// Create branches in the parameter tree and connect them to the variables above
	cout << endl << "Setting branch addresses for the new tree." << endl;
	param_tree->Branch("ped_near", &ped_near, "ped_near/D");
	param_tree->Branch("ped_far", &ped_far, "ped_far/D");
	param_tree->Branch("amp_near", &amp_near, "amp_near/D");
	param_tree->Branch("amp_far", &amp_far, "amp_far/D");
	param_tree->Branch("t_max_near", &t_max_near, "t_max_near/D");
	param_tree->Branch("t_max_far", &t_max_far, "t_max_far/D");
	param_tree->Branch("t_50_near", &t_50_near, "t_50_near/D");
	param_tree->Branch("t_50_far", &t_50_far, "t_50_far/D");
	
	// Variables for fitting 
	int n_bad_data=0,
	    n_bad_fit=0,
		n_fit=0;
	
	// ROOT histograms that will be filled with data from each event
	// The first stores the difference between the near and far signals' 50% rising edge times
	// The second holds the difference between the two signals' maximium times
	TH1D* h_t_diff_50=new TH1D("h_t_diff_50", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);
	TH1D* h_t_diff_max=new TH1D("h_t_diff_max", "Distribution of Time Difference; #Delta t  [ns]; count", 100, -50., 150.);
	
	// Iterates through each event in the data file. For take103, there are 3000 events
	cout << "Starting loop" << endl;
	for (int i=0; i<n_tot; i++){
		// Read from the raw data 
		raw_tree->GetEntry(i);
		
		// Find the normalization parameters and modify the variable values.
		// Note how the parameters are passed by reference
		int flag_near=channel_norm_fit(adc_near, &ped_near, &amp_near, &t_max_near, &t_50_near);
		int flag_far=channel_norm_fit(adc_far, &ped_far, &amp_far, &t_max_far, &t_50_far);
		cout << flag_near << " " << flag_far<<endl;
		
		// Transfer the parameter variables values into the branches of the appropriate tree
		param_tree->Fill();
		
		// Store additional data in the good events (those with no empty channels, etc.) 
		cout << is_good_far << endl;
		if (!is_good_near||!is_good_far){
			n_bad_data++;
		}
		else{
			// Add the timing details to the histograms only if the fit was acceptable
			// Otherwise, do nothing with the 50% rising edge times and maximum times
			if (flag_near==0&&flag_far==0){
				n_fit++;
				h_t_diff_50->Fill(t_50_far-t_50_near);
				h_t_diff_max->Fill(t_max_far-t_max_near);
			}
			
			else{
				n_bad_fit++;
			}
		}
		
		// Updates to the console every 500 events
		if ((i+1)%500==0) cout<<"\t\t"<<(i+1)<<" events have been processed."<<endl;
	}
	
	// Updates to the console
	cout<<"    "<<n_bad_data<<" events have empty channel(s)."<<endl;
	cout<<"    "<<n_bad_fit<<" events have bad fitting in either channel."<<endl;
	cout<<"    "<<n_fit<<" events have proper fits and time difference of signals were calculated."<<endl;
	
	// Change the directory to the parameter file we created and write the parameter tree to it
	param_file->cd();
	param_tree->Write();
	
	// Create plots for the two time-difference histograms.
	string plot_name = filePathPrefix  + fileName + "_normalization_parameters.pdf";
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
