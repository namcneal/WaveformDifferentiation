#include "../include/align_waveforms.h"

using namespace std;

/* Takes the original waveform data and applies the normalization parameters to the "good data"
 * Uses the raw data and normalization parameter ROOT files
 * 
 * Returns a new 'aligned' waveform. It has been aligned so that the waveform array starts a fixed distance
 * before the 50% rising edge. It has also been aligned vertically by subtracting the pedestal.
 */
 
int main(int argc, char **argv){
	string filePathPrefix = argv[1],
		   fileName = argv[2];
	
	string param_file_name = filePathPrefix + '/' + fileName + "_normalization_params.root",
		   param_tree_name = fileName + "_normalization_params",
		   raw_file_name = filePathPrefix + '/' + fileName + "_raw_data.root",
		   raw_data_tree_name = fileName + "_raw";
	 
	TFile* raw_file   = new TFile(raw_file_name.c_str(), "READ");
	TTree* raw_tree   = (TTree*)raw_file->Get(raw_data_tree_name.c_str());
  
	TFile* param_file = new TFile(param_file_name.c_str(), "READ");
	TTree* param_tree = (TTree*)param_file->Get(param_tree_name.c_str());

	if (!raw_tree){
		cout << "\n\nCannot find raw data tree" << endl << endl;
	}
	
	cout << endl << endl << param_file_name << endl << endl;
	if (!param_tree){
		cout << "Could not creat the new " << fileName << "_normalization_params tree." << endl;
	}
	
	double ped_near = 0, 
		   ped_far = 0, 
		   amp_near = 0, 
		   amp_far = 0, 
		   t_max_near = 0, 
		   t_max_far = 0, 
		   t_50_near = 0, 
		   t_50_far = 0;
	
	bool is_good_near = false, 
		 is_good_far  = false;
		
	int adc_near[N_width*32] = {0},
		adc_far[N_width*32] = {0};
		
	cout << "Setting branch addresses for the raw tree" << "...adc n";
	raw_tree->SetBranchAddress("adc_near", adc_near);
	cout << "...adc f";
	raw_tree->SetBranchAddress("adc_far", adc_far);
	cout << "...n quality";
	raw_tree->SetBranchAddress("near_quality", &is_good_near);
	cout << "...f quality." << endl;
	raw_tree->SetBranchAddress("far_quality", &is_good_far);
	
	cout << endl << "Setting branch addresses for the new tree." << endl;
	param_tree->SetBranchAddress("ped_near",   &ped_near);
	param_tree->SetBranchAddress("ped_far",    &ped_far);
	param_tree->SetBranchAddress("amp_near",   &amp_near);
	param_tree->SetBranchAddress("amp_far",    &amp_far);
	param_tree->SetBranchAddress("t_max_near", &t_max_near);
	param_tree->SetBranchAddress("t_max_far",  &t_max_far);
	param_tree->SetBranchAddress("t_50_near",  &t_50_near);
	param_tree->SetBranchAddress("t_50_far",   &t_50_far);
	
	cout << "Geting the total number of entries." << endl;
	int n_tot = raw_tree->GetEntries();
	cout << "Read "<< n_tot <<" events from " << fileName << "_raw_data.root" << endl;
	
	
	int aligned_adc_far[N_width*32] = {0};
	bool isNeutron = false;
	
	string aligned_waveform_file_name = filePathPrefix + '/' + fileName + "_aligned_waveforms.root",
		   aligned_data_tree_name = fileName + "_aligned_waveforms";
	TFile* aligned_file = new TFile(aligned_waveform_file_name.c_str(), "recreate");
	TTree* aligned_tree = new TTree(aligned_data_tree_name.c_str(), (fileName+"_aligned_waveforms").c_str());
	aligned_tree->Branch("aligned_ADC_far", aligned_adc_far, Form("adc_near[%i]/I", N_width*32));
	aligned_tree->Branch("isNeutron", &isNeutron, "isNeuton/O");
	
	int num_neutrons = 0,
		num_photons  = 0;
	cout << "Starting loop" << endl;
	for (int i=0; i<n_tot; i++){
		param_tree->GetEntry(i);
		raw_tree->GetEntry(i);
		
		if (is_good_far && is_good_near && t_max_far!=0. && t_max_near!=0.){
			align_waveform(adc_far, 
						   t_50_near, t_50_far, ped_far,
						   aligned_adc_far, isNeutron,
						   num_neutrons, num_photons);
						   
			aligned_tree->Fill();
		}
		
		if ((i+1)%500==0){
			cout<<"        "<<(i+1)<<" events have been processed."<<endl;
		}
	}
	
	cout << endl << "\t" << num_neutrons << " neutron events are present."<< endl
				 << "\t" << num_photons  << " photon events are present." << endl;
	
	aligned_file->cd();
	aligned_tree->Write();
	aligned_file->Close();
	return 0;
}


void align_waveform(int raw_adc_far[],
					int t_50_near, int t_50_far, int pedestal_far,
					int aligned_adc_far[], bool &isNeutron,
					int &num_neutrons, int &num_photons){
	
		double start_time = t_50_far - 20.;
		
		int start_TDC = (int)ceil(start_time);
		if (start_TDC%2){
			start_TDC += 1;
		}
		start_TDC/=2;
		
		// Fill in a new waveform array that starts at the index we chose with the 50% rising edge
		// The rest of the new array will be zeros
		for(int k = start_TDC; k < 32*N_width - start_TDC; k++){
			// Set the new array using the old array
			aligned_adc_far[k - start_TDC] = raw_adc_far[k];
			
			// Pedestal/background subtract the new array to align vertically
			aligned_adc_far[k - start_TDC] -= pedestal_far;
		}
		
		
		double t_50_diff = t_50_far - t_50_near;
		if (t_50_diff > neutron_low && t_50_diff < neutron_high){
			isNeutron = true;
			num_neutrons++;
		}
		else if (t_50_diff > photon_low && t_50_diff < photon_high){
			isNeutron = false;
			num_photons++;
		}
}


