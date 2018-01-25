#include "../include/generate_avg_waveform.h"
using namespace std;
					  
/* Takes the raw data and normalization parameter ROOT files.
 * Iterates through each event and adds each waveform to one of two 2-D histograms, the neutron or the gamma.
 * Gets the average neutron and average gamma waveforms from the X-profiles of the 2-D histograms. 
 * 
 * Arguments: dataFolderPath -> The path to the raw data file (example: "../../data") and 
     	      dataFileName   -> the name of the data file itself, without the extension (example: "take103")
 */
					  
int generate_avg_waveform(string dataFolderPath, string dataFileName){
	int numBins = 1000;
	
	cout << endl << endl << "Starting Program" << endl;
	generateAverageWaveforms(dataFolderPath, dataFileName, numBins);

	return 0;
}

void generateAverageWaveforms(string dataFolderPath, string dataFileName, int numBins){
	cout << "Reading in the data and creating the trees" << endl;
	
	// The names for reading form the file containing the normalization parameters
	string normalized_file_name = dataFolderPath + '/' + dataFileName + "_normalized_waveforms.root",
	       normalized_tree_name = dataFileName + "_normalized_waveforms";
  
	// Create the new file for the parameters and the 
	TFile* normalized_file = new TFile(normalized_file_name.c_str(), "READ");
	TTree* normalized_tree = (TTree*)normalized_file->Get(normalized_tree_name .c_str());
	
	// Check to make sure the proper files and trees are created. If not then a segmentation violation will occur. 
	if (!normalized_file) cout << endl << "Cannot find normalized waveform data file" << endl;
	if (!normalized_tree) cout << endl << "Cannot find normalized waveform data tree" << endl;
	
	cout << "Geting the total number of entries." << endl;
	int n_tot = normalized_tree->GetEntries();
	cout << "Read "<< n_tot <<" events from " << dataFileName << "_normalized_waveforms.root" << endl;
	
	// Variables to hold the aligned data
	int normalized_adc_event[N_width*32] = {0};  // Array to hold an aligned event
	bool isNeutron = false;						 // Whether the event was a neutron (true) or a gamma event (false)
	double t_50_rising_index = -1;			     // Time of the 50% rising edge, given as an index
	int start_index = -1;

	// Create branches in the aligned data tree and connect them to the variables above
	cout << endl << "Setting branch addresses for the normalized waveform data tree." << endl;
	normalized_tree->SetBranchAddress("normalized_adc_event", normalized_adc_event);
	normalized_tree->SetBranchAddress("isNeutron", &isNeutron);
	normalized_tree->SetBranchAddress("t_50_rising_index",  &t_50_rising_index);
	normalized_tree->SetBranchAddress("start_index",  &start_index);
	

	// 2D Histograms what will store 
	TH2D* hist_n=new TH2D("hist_n","Normalized Waveform: Neutron; time [ns]; Percentage",numBins,-20., 200.,140,-.2,1.2);
	TH2D* hist_p=new TH2D("hist_p","Normalized Waveform: Photon; time [ns]; Percentage",numBins,-20., 200.,140,-.2,1.2);
	
	int n_neutron=0;
	int n_photon=0;
	
	for (int i=0; i<n_tot; i++){
		normalized_tree->GetEntry(i);
		
		// Align the event according to the 50% rising edge and store it in the proper histogram
		if(isNeutron){
			for (int j=0; j<110; j++){
				hist_n->Fill(2.*(double)(start_index + j) - t_50_rising_index, normalized_adc_event[start_index + j]);
			}
			n_neutron++;
		}
		else{
			for (int j=0; j<110; j++){
				hist_p->Fill(2.*(double)(start_index + j) - t_50_rising_index , normalized_adc_event[start_index]);
			}
			n_photon++;
		}

		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"\t"<<n_neutron<<" neutron events are present."<<endl;
	cout<<"\t"<<n_photon<<" photon events are present."<<endl;
	
	// Create a file to store the averaged waveforms
	TFile* waveformFile = new TFile((dataFolderPath +  "/averaged_waveforms/" +  dataFileName + "_average_waveforms.root").c_str(), "RECREATE");

	// Write the original histograms to the file
	hist_n->Write();
	hist_p->Write();

	// Calculate the averages of the two histograms to obtain the averaged waveform and store the average to the file
	hist_n->ProfileX()->Write();
	hist_p->ProfileX()->Write();

	// Close all files
	waveformFile->Close();
	normalized_file->Close();
}
