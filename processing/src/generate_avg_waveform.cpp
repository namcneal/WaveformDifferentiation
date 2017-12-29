#include "generate_avg_waveform.h"
using namespace std;
					  
/* Takes the raw data and normalization parameter ROOT files.
 * Iterates through each event and adds each waveform to one of two 2-D histograms, the neutron or the gamma.
 * Gets the average neutron and average gamma waveforms from the X-profiles of the 2-D histograms. 
 */
					  
int main(int argc, char **argv){
	string dataFileName = argv[2],
		   dataFolderPath = argv[1];

	int numBins = 1000;
	generateAverageWaveforms(dataFolderPath, dataFileName, numBins);

	return 0;
}

void generateAverageWaveforms(string dataFolderPath, string dataFileName, int numBins){
	cout << "Reading in the data and creating the trees" << endl;
			   
	string plot_name = dataFolderPath + dataFileName + "PLEASE.pdf";
	
	TFile* norm_params_file=new TFile((dataFolderPath + dataFileName + "_normalization_params.root").c_str(), "READ");
	TTree* norm_params_tree=(TTree*)norm_params_file->Get((dataFileName + "_normalization_params").c_str());
	
	TFile* raw_file=new TFile((dataFolderPath + dataFileName+"_raw_data.root").c_str(), "READ");
	TTree* raw_tree=(TTree*)raw_file->Get((dataFileName + "_raw").c_str());
	
	double ped_n, ped_f, amp_n, amp_f, t_max_n, t_max_f, t_50_n, t_50_f;
	int is_good_n, is_good_f;
	int adc_n[N_width*32];
	int adc_f[N_width*32];
	
	cout << "Setting branch addresses to read into." << endl;
	norm_params_tree->SetBranchAddress("ped_n", &ped_n);
	norm_params_tree->SetBranchAddress("ped_f", &ped_f);
	norm_params_tree->SetBranchAddress("amp_n", &amp_n);
	norm_params_tree->SetBranchAddress("amp_f", &amp_f);
	norm_params_tree->SetBranchAddress("t_max_n", &t_max_n);
	norm_params_tree->SetBranchAddress("t_max_f", &t_max_f);
	norm_params_tree->SetBranchAddress("t_50_n", &t_50_n);
	norm_params_tree->SetBranchAddress("t_50_f", &t_50_f);
	
	raw_tree->SetBranchAddress("adc_n", adc_n);
	raw_tree->SetBranchAddress("adc_f", adc_f);
	raw_tree->SetBranchAddress("n_quality", &is_good_n);
	raw_tree->SetBranchAddress("f_quality", &is_good_f);
	
	int n_tot=raw_tree->GetEntries();
	cout<<"Read "<<n_tot<<" events from "<<dataFileName<<"_raw_data.root"<<endl;
	
	TH2D* hist_n=new TH2D("hist_n","Normalized Waveform: Neutron; time [ns]; Percentage",numBins,-20., 200.,140,-.2,1.2);
	TH2D* hist_p=new TH2D("hist_p","Normalized Waveform: Photon; time [ns]; Percentage",numBins,-20., 200.,140,-.2,1.2);
	
	int n_neutron=0;
	int n_photon=0;
	
	for (int i=0; i<n_tot; i++){
		norm_params_tree->GetEntry(i);
		raw_tree->GetEntry(i);
		
		if (is_good_f==1 && is_good_n==1 && t_max_f!=0. && t_max_n!=0.){
			double t_50_diff=t_50_f-t_50_n;
			double start_time=t_50_f-20.;
			int start_TDC=(int)ceil(start_time);
			if (start_TDC%2) start_TDC+=1;
			start_TDC/=2;
			
			if (t_50_diff>neutron_low && t_50_diff<neutron_high){
				for (int j=0; j<110; j++){
					hist_n->Fill(2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
				}
				n_neutron++;
			}
			if (t_50_diff>photon_low && t_50_diff<photon_high){
				for (int j=0; j<110; j++){
					hist_p->Fill(2.*(double)(start_TDC+j)-t_50_f,(adc_f[start_TDC+j]-ped_f)/amp_f);
				}
				n_photon++;
			}
		}
		if ((i+1)%500==0) cout<<"        "<<(i+1)<<" events have been processed."<<endl;
	}
	cout<<"    "<<n_neutron<<" neutron events are present."<<endl;
	cout<<"    "<<n_photon<<" photon events are present."<<endl;
	
	TFile* waveformFile = new TFile((dataFolderPath + dataFileName + "_average_waveforms.root").c_str(), "RECREATE");
	hist_n->ProfileX()->Write();
	hist_p->ProfileX()->Write();
	waveformFile->Close();
	
	norm_params_file->Close();
	raw_file->Close();
}