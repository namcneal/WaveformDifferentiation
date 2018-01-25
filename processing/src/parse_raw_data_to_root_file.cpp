#include "../include/parse_raw_data_to_root_file.h"
using namespace std;

/* Takes a raw .dat data file and tranfers the contents into a .root file of the same name.
 * 
 * Use as parser <filename(-".dat")> <N_width> <N_spill> <channel_1> "n/f" <channel_2> "n/f"
 * Requires channel_1 < channel_2
 * 
 * Example for additional arguments after data file name: 15 200 0 n 3 f 
 * 
 * The N_width, N_spill, and other arguments are file specific and depend on the run. 
*/
int main(int argc, char **argv){
	// The path to the raw data file (example: ../../data) and 
     	// the name of the data file itself, without the extension (example: take103)
	string filePathPrefix = argv[1],
		   fileName   = argv[2];
	
	// Parameters for the parsing 
	int n_width=atoi(argv[3]);	// Number of words in the ADC event 
	const int sampling_len=n_width*32;	// Total length of each ADC event 
	int n_spill=atoi(argv[4]); 	// Number of spills
	
	int chan1=atoi(argv[5]);	// 
	char* label1=argv[6];		// Label for channel 1, near or far
	int chan2=atoi(argv[7]);	// 
	char* label2=argv[8];		// Label for channel 2, near or far
	
	int n_event=2048*16/(sampling_len*4+HEADER_LEN);
	int n_skip=2048*16-n_event*(sampling_len*4+HEADER_LEN);
	int index=0;
	int junk, buffer;
	char junk_string[20];
	
	// Define the full path to the data file using the passed parameters/arguments
	string file_name = filePathPrefix + '/' + fileName + ".dat";
	cout << file_name << endl << endl;
	
	// Open the data file 
	FILE* in;
	in = fopen(file_name.c_str(), "r");
	if(!in){
		cout << endl << "Could not open data (.dat) file" << endl;
	}
	
	// Create the ROOT file that will store the raw data that we transfer from the .dat file
	string root_name = filePathPrefix + '/' + fileName + "_raw_data.root";
	TFile *file = new TFile (root_name.c_str(), "recreate");

	// Create a tree inside the ROOT file
	string tree_name = fileName + "_raw";
	TTree *tree =  new TTree (tree_name.c_str(), tree_name.c_str());
	
	// Variables for reading in data
	int adc_near[sampling_len]; 	// Stores the events from the near detector
	int adc_far[sampling_len];	    // Stores the events from the far detector
	
	for(int i = 0; i < sampling_len; i++){
		adc_near[i] = 0;
		adc_far[i]  = 0;
	}
	
	// Variables for judging the quality of an event signal
	bool near_quality = 0,
	     far_quality = 0;
	
	int num_near_good = 0,
	    num_far_good = 0;

	// Create branches in the raw data tree to store data for each event.
	// Four branches for each of the variables defined above. 	
	tree->Branch("adc_near", adc_near, Form("adc_near[%i]/I", sampling_len));
	tree->Branch("adc_far", adc_far, Form("adc_far[%i]/I", sampling_len));
	tree->Branch("near_quality", &near_quality, "near_quality/O");
	tree->Branch("far_quality", &far_quality, "far_quality/O");
	
	
	// Iterate through each event in the data file and calculate each of the four variables 
	for (int i=0; i<n_spill; i++){
		for (int j=0; j < n_event; j++){
			
			// Read in the header and then dump it
			for (int k=0; k<HEADER_LEN; k++) {
				fscanf(in, "%x", &junk);
			}

			// Read through the relevant contents of each event
			for (int k=0; k<sampling_len; k++){

				// Read in data from Channel 1
				for (int l=0; l<chan1; l++){
					fscanf(in, "%x", &junk);
				}
				
				// Write/transfer the Channel 1 data into either the near or the far array
				if ((strcmp(label1,"n")==0) && (strcmp(label2,"f")==0)){
					fscanf(in, "%x", &buffer);
					adc_near[k]=(buffer-32768);
				}
				else if ((strcmp(label1,"f")==0) && (strcmp(label2,"n")==0)){
					fscanf(in, "%x", &buffer);
					adc_far[k]=(buffer-32768);
				}
				
				// Read in the data from Channel 2
				for (int l=0; l<(chan2-chan1-1); l++){
					fscanf(in, "%x", &junk);
				}
				
				// Write/transfer the Channel 2 data into either the near or the far array
				if ((strcmp(label1,"n")==0)&&(strcmp(label2,"f")==0)){
					fscanf(in, "%x", &buffer);
					adc_far[k]=(buffer-32768);
				}
				else if ((strcmp(label1,"f")==0)&&(strcmp(label2,"n")==0)){
					fscanf(in, "%x", &buffer);
					adc_near[k]=(buffer-32768);
				}

				for (int l=0; l<(3-chan2); l++){fscanf(in, "%x", &junk);}
			}

			// Excludes j = 0 to get rid off the mis-triggered first event in each spill
			if (j != 0) {
				// Check the near and far events
				near_quality = quality_check(adc_near, sampling_len);
				far_quality = quality_check(adc_far, sampling_len);

				// Make a note of the number of good events for both near and far
				if (near_quality) num_near_good++;
				if (far_quality)  num_far_good++;
				
				// Fill the next spot in the ROOT tree with the values we have caluclated for the four variables
				tree->Fill();

				// Next event
				index++;
			}
		}
		
		// Updates
		for (int l=0; l<n_skip; l++){fscanf(in, "%x", &junk);}
		fscanf(in, "%i", &junk);
		for (int l=0; l<7; l++){fscanf(in, "%s", junk_string);}
	}

	// Updates
	cout<<index<<" events are parsed out of "<< n_spill<<" spills."<<endl;
	cout<<"    "<< num_near_good << " events have signal in the near channel."<<endl;
	cout<<"    "<< num_far_good  << " events have signal in the far channel."<<endl;
	tree->Write();
	fclose(in);

	file->Close();
	return 0;
}

// Takes an event array and its length
// Returns true  for good data and 
//         false for empty channel
bool quality_check(int* adc, int sampling_len){
    bool flag = false;
    
	for (int i = 0; i < sampling_len; i++){
		if (adc[i] > 400) {
			flag = true; 
			break;
		}
	}
	
	return flag;
}
