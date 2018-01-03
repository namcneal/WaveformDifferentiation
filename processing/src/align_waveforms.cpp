#include "../include/align_waveforms.h"

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

	// Create a new ROOT file with a tree for the aligned events
	string aligned_waveform_file_name = filePathPrefix + '/' + fileName + "_aligned_waveforms.root",
	       aligned_data_tree_name = fileName + "_aligned_waveforms";
	TFile* aligned_file = new TFile(aligned_waveform_file_name.c_str(), "recreate");
	TTree* aligned_tree = new TTree(aligned_data_tree_name.c_str(), (fileName+"_aligned_waveforms").c_str());

	// Define variables for the aligned data
	int adc_event_data[N_width*32] = {0};   // Array to hold an aligned event
	bool isNeutron = false;			// Whether the event was a neutron (true) or a gamma event (false)
	int t_50_rising_edge = -1;		// Time of the 50% rising edge, given as an index
	
	// Create branches in the tree to store the data for each event
	aligned_tree->Branch("ADC_event_data", adc_event_data, Form("adc_event_data[%i]/I", N_width*32));
	aligned_tree->Branch("isNeutron", &isNeutron, "isNeuton/O");
	aligned_tree->Branch("t_50_rising_edge",  &t_50_rising_edge,  "t_50_rising_edge");
	
	// Variables for keeping track, not stored 
	int num_neutrons = 0,
	    num_photons  = 0;
	
	// Iterate through each event to subtract the pedestal and shift it based on 50% rising edge
	cout << "Starting loop" << endl;
	for (int i=0; i<n_tot; i++){

		// Read from the two trees
		raw_tree->GetEntry(i);
		param_tree->GetEntry(i);

		// Align and store the waveforms only for the good events for which the maximum time is positive as well
		if (is_good_far && is_good_near && t_max_far!=0. && t_max_near!=0.){
			// Perform the vertical (pedestal correction) and horizontal (time shifting) alignment			
			align_waveform(adc_far, 
					   t_50_near, t_50_far, ped_far, amp_far,
					   aligned_adc_far, isNeutron,
					   num_neutrons, num_photons);

			// Transfer the 50% edge index to the new file
			t_50_rising_edge = t_50_far;

			// Transfer the variable values into the tree		   
			aligned_tree->Fill();
		}
		
		// Updates
		if ((i+1)%500==0) cout << "\n\n" << (i+1) << " events have been processed."<<endl;
		
	}
	
	cout << endl << "\t" << num_neutrons << " neutron events are present."<< endl
				 << "\t" << num_photons  << " photon events are present." << endl;
	
	// Write the tree with the aligned events to the proper ROOT file
	aligned_file->cd();
	aligned_tree->Write();
	aligned_file->Close();
	return 0;
}


void align_waveform(int raw_adc_far[],
			int t_50_near, int t_50_far, int pedestal_far, double amp_far, 
			int aligned_adc_far[], bool &isNeutron,
			int &num_neutrons, int &num_photons){
	
	// Background/pedestal subtraction and then normalize
	for(int k = start_TDC; k < 32*N_width; k++){
		// Set the new array using the old array
		aligned_adc_far[k] = raw_adc_far[k];
	
		// Pedestal/background subtract the new array to align vertically
		aligned_adc_far[k] -= pedestal_far;

		// Normalize based on amplification parameter
		aligned_adc_far[k] /= amp_far
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


