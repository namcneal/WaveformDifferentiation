#include "../qval_calculation_functions.h"

string dataFolderPath = "../data",
       dataFileName = "take103";

string avgWaveformROOTFileName = dataFolderPath + "/" + dataFileName + "_average_waveforms.root";

TFile * waveformFile = new TFile(avgWaveformROOTFileName.c_str(), "READ");
    
if(!waveformFile){
    cout << "Could not find/read the waveform file";
}

cout << "Read in the waveform file";

TProfile * avgNeutronWaveform = NULL,
         * avgPhotonWaveform  = NULL;
waveformFile->GetObject("hist_n_pfx", avgNeutronWaveform);
waveformFile->GetObject("hist_p_pfx", avgPhotonWaveform);

#define numSteps 150

int minTailStart = 10,
    maxTailStart  = 100,
    tailEnd      = 170;

double stepSize = (maxTailStart - minTailStart) / numSteps;

double tailStart = 0,
       neutronQVal = 0,
       photonQVal  = 0;

double tailStartArray[numSteps],
       neutronQValArray[numSteps],
       photonQValArray[numSteps];

int i = 0;
for(tailStart = minTailStart; tailStart <= MaxTailStart; tailStart += stepSize){
    
    cout << i << endl;
    calculateQVals(tailStart, tailEnd, avgNeutronWaveform, avgPhotonWaveform,
                   neutronQVal, photonQVal);
        
    tailStartArray[i] = tailStart;
    neutronQValArray[i] = neutronQVal;
    photonQValArray[i]  = photonQVal;
    
    i++;
}