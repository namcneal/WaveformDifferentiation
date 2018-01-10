#include <iostream>
#include <string>

#include "TFile.h"
#include "TProfile.h"

using namespace std;

void calculateQVal(double tailStart, double tailEnd, 
				   TProfile* neutronWaveform, TProfile* photonWaveform,
				   double &neutronQVal, double &photonQVal);

int main(int argc, char** argv){
	
	string dataFolderPath = argv[1],
		   dataFileName   = argv[2];
	   
	string avgWaveformROOTFileName = dataFolderPath + "/" + dataFileName + "_average_waveforms.root";

	TFile * waveformFile = new TFile(avgWaveformROOTFileName.c_str(), "READ");

	if(!waveformFile){
		cout << "Could not find/read the waveform file" << endl;
	}

	cout << "Read in the waveform file";

	TProfile * avgNeutronWaveform = NULL,
			 * avgPhotonWaveform  = NULL;
	waveformFile->GetObject("hist_n_pfx", avgNeutronWaveform);
	waveformFile->GetObject("hist_p_pfx", avgPhotonWaveform);

	int numStepsStart = 25,
		minTailStart  = 30,
		maxTailStart  = 100,
		
		numStepsEnd   = 25,
		minTailEnd    = 120,
		maxTailEnd    = 200;
		
	double stepSizeStart = (maxTailStart - minTailStart) / numStepsStart,
		   stepSizeEnd   = (maxTailEnd   - minTailEnd  ) / numStepsEnd;

	double tailStart = 0,
		   tailEnd   = 0,
		   neutronQVal = 0,
		   photonQVal  = 0;

//	double tailStartArray[numSteps],
//		   tailEndArray[numSteps];

	double neutronQVals[numStepsStart][numStepsEnd],
		   photonQVals[numStepsStart][numStepsEnd];

	int i = 0,
		j = 0;
	for(tailStart = minTailStart; tailStart <= maxTailStart; tailStart += stepSizeStart){
		for(tailEnd = minTailEnd; tailEnd <= maxTailEnd; tailEnd += stepSizeEnd){
		
			calculateQVal(tailStart, tailEnd, 
						  avgNeutronWaveform, avgPhotonWaveform,
						  neutronQVal, photonQVal);
				
			neutronQVals[i][j] = neutronQVal;
			photonQVals[i][j]  = photonQVal;
		
			j++;
		}
		if(i%20 ==0) cout << i << endl;
		i++;
	}
		
	
}

void calculateQVal(double tailStart, double tailEnd, 
				   TProfile* avgNeutron, TProfile* avgPhoton,
				   double &neutronQVal, double &photonQVal){
	
	int tailStartBin = avgNeutron->FindBin(tailStart),
		tailEndBin   = avgNeutron->FindBin(tailEnd);
		
	double a_main_n=0.,
		   a_tail_n=0.,
		   a_main_p=0.,
		   a_tail_p=0.;
	
	for (int i=1; i<tailStartBin; i++){
		a_main_n+=avgNeutron->GetBinContent(i);
		a_main_p+=avgPhoton->GetBinContent(i);
	}
	for (int i = tailEndBin; i<=avgNeutron->FindBin(199.75); i++){
		a_tail_n+=avgNeutron->GetBinContent(i);
		a_tail_p+=avgPhoton->GetBinContent(i);
	}
	
	double neutronCurrentQVal = a_tail_n / (a_tail_n+a_main_n),
		   photonCurrentQVal  = a_tail_p / (a_tail_p+a_main_p);
		   
	neutronQVal = neutronCurrentQVal;
	photonQVal  = photonCurrentQVal;
}
