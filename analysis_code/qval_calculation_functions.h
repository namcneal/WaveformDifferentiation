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
# include "TH1.h"
# include "TH2.h"
# include "TH1F.h"
# include "TF1.h"
# include "TApplication.h"
# include "TMultiGraph.h"
# include "TMath.h"
# include "TTree.h"
# include "TFile.h"
# include "TLine.h"
# include "TNtuple.h"
# include "TProfile.h"

void calculateQVal(int tailStart, int tailEnd, 
				   TProfile* neutronWaveform, TProfile* photonWaveform,
				   double &neutronQVal, double &photonQVal);
				   