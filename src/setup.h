#ifndef SETUP_H
#define SETUP_H

const Int_t nRealCuts = 6;//Number of cuts made
const Int_t nMCCuts = 5;//Number of cuts made

void HistArraySetupReal(TH1D **hist, TH1D *hCut[][nRealCuts], Int_t nBins,
		   Double_t min,
		   Double_t max, Int_t iCut, TString hName);


void FillCutsReal(TH1D *hist[][nRealCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist);

void HistArraySetupMC(TH1D **hist, TH1D *hCut[][nMCCuts], Int_t nBins,
		   Double_t min,
		   Double_t max, Int_t iCut, TString hName);


void FillCutsMC(TH1D *hist[][nMCCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist);

#endif
