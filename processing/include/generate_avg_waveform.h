# include <iostream>
# include <string>

# include "TROOT.h"
# include "TTree.h"
# include "TFile.h"
# include "TProfile.h"
# include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"


const double PERTDC=2.;
const int N_width=15;
const double photon_low=-20.;
const double neutron_low=30.;
const double photon_high=20.;
const double neutron_high=150.;

void generateAverageWaveforms(std::string dataFolderPath, std::string dataFileName, int numBins);
