#include "../include/normalize_waveforms.h"

using namespace std;

/* Takes the original waveform data and applies the normalization parameters to the "good data" of the FAR channel. 
 * Uses the raw data and normalization parameter ROOT files
 * 
 * Returns a ROOT file with the new normalized waveforms. It has been baseline subtracted and then normalized through the amplification parameter. 
 */
 
int main(int argc, char **argv){
	
	// The path to the raw data file (example: ../../data) and 
     	// the name of the data file itself, without the extension (example: take103)
	string filePathPrefix = argv[1],
		   fileName = argv[2];
	
	// The names for reading in raw data from the ROOT file
	string raw_file_name = filePathPrefix + '/' + fileName + "_raw_data.root",
	       raw_data_tree_name = fileName + "_raw";
	
	// The names for reading form the file containing the normalization parameters
	string param_file_name = filePathPrefix + '/' + fileName + "_normalization_params.root",
	       param_tree_name = fileName + "_normalization_params";

	// Create the variables necessary for reading from the raw data ROOT file and its tree
	TFile* raw_file   = new TFile(raw_file_name.c_str(), "READ");
	TTree* raw_tree   = (TTree*)raw_file->Get(raw_data_tree_name.c_str());
  
	// Create the new file for the parameters and the 
	TFile* param_file = new TFile(param_file_name.c_str(), "READ");
	TTree* param_tree = (TTree*)param_file->Get(param_tree_name.c_str());
	
	// Check to make sure the proper files and trees are created. If not then a segmentation violation will occur. 
	if (!raw_file)   cout << endl << "Cannot find raw data file" << endl;
	if (!raw_tree)   cout << endl << "Cannot find raw data tree" << endl;
	if (!param_file) cout << endl << "Cannot find parameter data file" << endl;
	if (!param_tree) cout << endl << "Cannot find parameter data tree" << endl;
	
	// Create the variables that will be read into from the raw data ROOT file
	bool is_good_near = false, 
      	     is_good_far  = false;
		
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
	param_tree->SetBranchAddress("ped_near",   &ped_near);
	param_tree->SetBranchAddress("ped_far",    &ped_far);
	param_tree->SetBranchAddress("amp_near",   &amp_near);
	param_tree->SetBranchAddress("amp_far",    &amp_far);
	param_tree->SetBranchAddress("t_max_near", &t_max_near);
	param_tree->SetBranchAddress("t_max_far",  &t_max_far);
	param_tree->SetBranchAddress("t_50_near",  &t_50_near);
	param_tree->SetBranchAddress("t_50_far",   &t_50_far);

	// Create a new ROOT file with a tree for the normalized events
	string normalized_waveform_file_name = filePathPrefix + '/' + fileName + "_normalized_waveforms.root",
	       normalized_data_tree_name = fileName + "_normalized_waveforms";
	TFile* normalized_file = new TFile(normalized_waveform_file_name.c_str(), "recreate");
	TTree* normalized_tree = new TTree(normalized_data_tree_name.c_str(), (fileName+"_normalized_waveforms").c_str());

	// Define variables for the normalized data
	int normalized_adc_event[N_width*32] = {0};   // Array to hold an normalized event
	bool isNeutron = false;			// Whether the event was a neutron (true) or a gamma event (false)
	int t_50_rising_index = -1;		// Time of the 50% rising edge, given as an index
	int start_index = -1;			// Time of the start of an event, based on the 50% rising edge
	
	// Create branches in the tree to store the data for each event
	normalized_tree->Branch("normalized_adc_event", normalized_adc_event, Form("normalized_adc_event[%i]/I", N_width*32));
	normalized_tree->Branch("isNeutron", &isNeutron, "isNeuton/O");
	normalized_tree->Branch("t_50_rising_index",  &t_50_rising_index,  "t_50_rising_index");
	normalized_tree->Branch("start_index",  &start_index,  "start_index");
	
	
	// Variables for keeping track, not stored 
	int num_neutrons = 0,
	    num_photons  = 0;
	
	// Iterate through each event to subtract the pedestal and shift it based on 50% rising edge
	cout << "Starting loop" << endl;
	for (int i=0; i<n_tot; i++){

		// Read from the two trees
		raw_tree->GetEntry(i);
		param_tree->GetEntry(i);

		// normalize and store the waveforms only for the good events for which the maximum time is positive as well
		if (is_good_far && is_good_near && t_max_far!=0. && t_max_near!=0.){

			// Define the start of the event based on the 50% rising edge 
			double start_time = t_50_far - 20.;
			int start_TDC = (int)ceil(start_time);
			if (start_TDC%2){
				start_TDC += 1;
			}
			start_TDC/=2;

			// Set the starting index that will be transferred to the normalized waveform file
			start_index = start_TDC;

			// Perform the vertical (pedestal correction) and scaling based on amplitude		
			normalize_waveform(adc_far, 
					   t_50_near, t_50_far, ped_far, amp_far,
					   normalized_adc_event, isNeutron,
					   num_neutrons, num_photons);

			// Transfer the 50% edge index to the new file
			t_50_rising_index = t_50_far;

			// Transfer the variable values into the tree		   
			normalized_tree->Fill();
		}
		
		// Updates
		if ((i+1)%500==0) cout << "\n\n" << (i+1) << " events have been processed."<<endl;
		
	}
	
	cout << endl << "\t" << num_neutrons << " neutron events are present."<< endl
				 << "\t" << num_photons  << " photon events are present." << endl;
	
	// Write the tree with the normalized events to the proper ROOT file
	normalized_file->cd();
	normalized_tree->Write();
	normalized_file->Close();
	return 0;
}


void normalize_waveform(int raw_adc_far[],
			int t_50_near, int t_50_far, int pedestal_far, double amp_far, 
			int normalized_adc_far[], bool &isNeutron,
			int &num_neutrons, int &num_photons){
	
	// Background/pedestal subtraction and then normalize
	for(int k = 0; k < 32*N_width; k++){
		// Set the new array using the old array
		normalized_adc_far[k] = raw_adc_far[k];
	
		// Pedestal/background subtract the new array to align vertically
		normalized_adc_far[k] -= pedestal_far;

		// Normalize based on amplification parameter
		normalized_adc_far[k] /= amp_far;
	}

	// Tag the event as either a neutron or a gamma based on the time difference between near and far channels
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


