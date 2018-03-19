#ifndef SETUP_H
#define SETUP_H

const Int_t nRealCuts = 7;//Number of cuts made
const Int_t nMCCuts = 7;//Number of cuts made

//Real
///////////////
void HistArraySetupReal(TH1D **hist, TH1D *hCut[][nRealCuts], Int_t nBins,
		   Double_t min,
		   Double_t max, Int_t iCut, TString hName);

void Hist2D_ArraySetupReal(TH2D **hist, TH2D *hCut[][nMCCuts],
			 Int_t nBins_x, Double_t min_x, Double_t max_x,
			 Int_t nBins_y, Double_t min_y, Double_t max_y,
			 Int_t iCut, TString hName);

void FillCutsReal(TH1D *hist[][nRealCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist);

void Fill2D_CutsReal(TH2D *hist[][nRealCuts], Double_t *variables,
		     Int_t cut, Int_t nCutHist);


//MC
void HistArraySetupMC(TH1D **hist, TH1D *hCut[][nMCCuts], Int_t nBins,
		   Double_t min,
		   Double_t max, Int_t iCut, TString hName);

void Hist2D_ArraySetupMC(TH2D **hist, TH2D *hCut[][nMCCuts],
			 Int_t nBins_x, Double_t min_x, Double_t max_x,
			 Int_t nBins_y, Double_t min_y, Double_t max_y,
			 Int_t iCut, TString hName);

void FillCutsMC(TH1D *hist[][nMCCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist);

void Fill2D_CutsMC(TH2D *hist[][nMCCuts], Double_t *variables,
		     Int_t cut, Int_t nCutHist);

#endif
