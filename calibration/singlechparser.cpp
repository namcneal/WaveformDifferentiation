# include <cstdio>
# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
# include <cassert>
# include <vector>
# include <algorithm>

# include "TCanvas.h"
# include "TROOT.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH1F.h"
# include "TF1.h"
# include "TLegend.h"
# include "TLatex.h"
# include "TStyle.h"
# include "TApplication.h"
# include "TMultiGraph.h"
# include "TMath.h"
# include "TTree.h"
# include "TFile.h"
# include "TLine.h"
# include "TNtuple.h"

using namespace std;
int header_len=12;

///////////////////////////////////////////////////////
//
//Use as parser <filename(-".dat")> <N_width> <N_spill> <channel_1> 
//
///////////////////////////////////////////////////////

int quality_check(int* adc, int sampling_len){
	int flag=1;
	for (int i=0; i<sampling_len; i++){
		if (adc[i]>=4095) flag=0;
	}
	return flag;//flag==1 good data; flag==0 existing 'ffff' etc.
}

int main(int argc, char **argv){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gErrorIgnoreLevel = kWarning;
	
	string name (argv[1]);
	int n_width=atoi(argv[2]);
	int sampling_len=n_width*32;
	int n_spill=atoi(argv[3]);
	int chan1=atoi(argv[4]);
	
	int n_event=2048*16/(sampling_len*4+header_len);
	int n_skip=2048*16-n_event*(sampling_len*4+header_len);
	int index=0;
	int junk, buffer;
	char junk_string[20];
	
	string root_name=name+"_raw_data.root";
	FILE* in;
	in=fopen((name+".dat").c_str(), "r");
	TFile *file=new TFile (root_name.c_str(), "recreate");
	TTree *tree=new TTree (Form("%s_raw", name.c_str()),Form("%s_raw", name.c_str()));
	
	int adc_n[sampling_len];
	for (int i=0; i<sampling_len; i++){
		adc_n[i]=0;
	}
	int n_quality=0;
	int n_good=0;
	
	tree->Branch("adc_n", adc_n, Form("adc_n[%i]/I", sampling_len));
	tree->Branch("n_quality", &n_quality, "n_quality/I");
	
	for (int i=0; i<n_spill; i++){
		for (int j=0; j<n_event; j++){
			for (int k=0; k<header_len; k++) {fscanf(in, "%x", &junk);}
			for (int k=0; k<sampling_len; k++){
				for (int l=0; l<chan1; l++){fscanf(in, "%x", &junk);}
				fscanf(in, "%x", &buffer);
				adc_n[k]=(buffer-32768);
				for (int l=0; l<(3-chan1); l++){fscanf(in, "%x", &junk);}
			}
			if (j!=0) {//get rid off the mis-triggered first event in each spill
				n_quality=quality_check(adc_n, sampling_len);
				if (n_quality==1) n_good++;
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
	tree->Write();
	fclose(in);

	file->Close();
	return 0;
}
