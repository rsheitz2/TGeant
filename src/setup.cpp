#include "common.hxx"
#include "setup.h"

void HistArraySetupReal(TH1D **hist, TH1D *hCut[][nRealCuts], Int_t nBins,
			Double_t min,
			Double_t max, Int_t iCut, TString hName){
  
  for (Int_t i=0; i<nRealCuts; i++) {
    hist[i] = new TH1D(Form("hCut%i_%s", i, hName.Data() ),
		       Form("hCut%i_%s", i, hName.Data() ),
		       nBins, min, max);

    hCut[iCut][i] = hist[i];
  }
}//HistArraySetupReal


void FillCutsReal(TH1D *hist[][nRealCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist){
  for (Int_t i=0; i<nCutHist; i++) hist[i][cut]->Fill(variables[i] );
}//FillCutsReal


void HistArraySetupMC(TH1D **hist, TH1D *hCut[][nMCCuts], Int_t nBins,
			Double_t min,
			Double_t max, Int_t iCut, TString hName){
  
  for (Int_t i=0; i<nMCCuts; i++) {
    hist[i] = new TH1D(Form("hCut%i_%s", i, hName.Data() ),
		       Form("hCut%i_%s", i, hName.Data() ),
		       nBins, min, max);

    hCut[iCut][i] = hist[i];
  }
}//HistArraySetupMC


void FillCutsMC(TH1D *hist[][nMCCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist){
  for (Int_t i=0; i<nCutHist; i++) hist[i][cut]->Fill(variables[i] );
}//FillCutsMC
