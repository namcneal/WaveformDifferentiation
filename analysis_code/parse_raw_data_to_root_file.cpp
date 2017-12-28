#include "parser.h"
using namespace std;

/* Takes a raw .dat data file and tranfers the contents into a .root file of the same name.
 * 
 * Use as parser <filename(-".dat")> <N_width> <N_spill> <channel_1> "n/f" <channel_2> "n/f"
 * Requires channel_1 < channel_2
 * 
 * The N_width, N_spill, and other arguments are file specific and depend on the run. 
*/
int main(int argc, char **argv){
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	
	string filePathPrefix = argv[1],
		   fileName 	  = argv[2];
	
	
	int n_width=atoi(argv[3]);
	int sampling_len=n_width*32;
	int n_spill=atoi(argv[4]);
	int chan1=atoi(argv[5]);
	char* label1=argv[6];
	int chan2=atoi(argv[7]);
	char* label2=argv[8];
	int n_event=2048*16/(sampling_len*4+HEADER_LEN);
	int n_skip=2048*16-n_event*(sampling_len*4+HEADER_LEN);
	int index=0;
	int junk, buffer;
	char junk_string[20];
	
	string root_name = filePathPrefix + fileName + "_raw_data.root";
	FILE* in;
	in = fopen((filePathPrefix + fileName + ".dat").c_str(), "r");
	
	TFile *file = new TFile (root_name.c_str(), "recreate");
	
	string tree_name = fileName + "_raw";
	TTree *tree =  new TTree (tree_name.c_str(), tree_name.c_str());
	
	int adc_n[sampling_len];
	int adc_f[sampling_len];
	for (int i=0; i<sampling_len; i++){
		adc_n[i]=0;
		adc_f[i]=0;
	}
	int n_quality=0,
	    f_quality=0;
	
	int n_good=0;
	int f_good=0;
	
	tree->Branch("adc_n", adc_n, Form("adc_n[%i]/I", sampling_len));
	tree->Branch("adc_f", adc_f, Form("adc_f[%i]/I", sampling_len));
	tree->Branch("n_quality", &n_quality, "n_quality/I");
	tree->Branch("f_quality", &f_quality, "f_quality/I");
	
	
	for (int i=0; i<n_spill; i++){
		for (int j=0; j<n_event; j++){
			for (int k=0; k<HEADER_LEN; k++) {fscanf(in, "%x", &junk);}
			for (int k=0; k<sampling_len; k++){
				for (int l=0; l<chan1; l++){fscanf(in, "%x", &junk);}
				if ((strcmp(label1,"n")==0)&&(strcmp(label2,"f")==0)){
					fscanf(in, "%x", &buffer);
					adc_n[k]=(buffer-32768);
				}
				if ((strcmp(label1,"f")==0)&&(strcmp(label2,"n")==0)){
					fscanf(in, "%x", &buffer);
					adc_f[k]=(buffer-32768);
				}
				for (int l=0; l<(chan2-chan1-1); l++){fscanf(in, "%x", &junk);}
				if ((strcmp(label1,"n")==0)&&(strcmp(label2,"f")==0)){
					fscanf(in, "%x", &buffer);
					adc_f[k]=(buffer-32768);
				}
				if ((strcmp(label1,"f")==0)&&(strcmp(label2,"n")==0)){
					fscanf(in, "%x", &buffer);
					adc_n[k]=(buffer-32768);
				}
				for (int l=0; l<(3-chan2); l++){fscanf(in, "%x", &junk);}
			}
			if (j!=0) {//get rid off the mis-triggered first event in each spill
				n_quality = (int) quality_check(adc_n, sampling_len);
				if (n_quality == 1) n_good++;
				f_quality = (int) quality_check(adc_f, sampling_len);
				if (f_quality == 1) f_good++;
				tree->Fill();
				index++;
				}
		}
		for (int l=0; l<n_skip; l++){fscanf(in, "%x", &junk);}
		fscanf(in, "%i", &junk);
		for (int l=0; l<7; l++){fscanf(in, "%s", junk_string);}
	}
	cout<<index<<" events are parsed out of "<< n_spill<<" spills."<<endl;
	cout<<"    "<<n_good<<" events have signal in the near channel."<<endl;
	cout<<"    "<<f_good<<" events have signal in the far channel."<<endl;
	tree->Write();
	fclose(in);

	file->Close();
	return 0;
}

bool quality_check(int* adc, int sampling_len){//return 1, good data; return 0, empty channel
    bool flag = false;
    
	for (int i=0; i<sampling_len; i++){
		if (adc[i]>400) {flag = true; break;}
	}
	
	return flag;
}
