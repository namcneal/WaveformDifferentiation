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

const int HEADER_LEN = 12;

/* Used to check the ADC data for empty channels.
 * Returns true if the data is good and
 * 		   false for empty channels
*/
bool quality_check(int* adc, int sampling_len);