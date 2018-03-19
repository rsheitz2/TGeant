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


void Hist2D_ArraySetupReal(TH2D **hist, TH2D *hCut[][nRealCuts],
			   Int_t nBins_x, Double_t min_x, Double_t max_x,
			   Int_t nBins_y, Double_t min_y, Double_t max_y,
			   Int_t iCut, TString hName){
  
  for (Int_t i=0; i<nRealCuts; i++) {
    hist[i] = new TH2D(Form("hCut%i_%s", i, hName.Data() ),
		       Form("hCut%i_%s", i, hName.Data() ),
		       nBins_x, min_x, max_x,
		       nBins_y, min_y, max_y);

    hCut[iCut][i] = hist[i];
  }
}//Hist2D_ArraySetupReal


void FillCutsReal(TH1D *hist[][nRealCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist){
  for (Int_t i=0; i<nCutHist; i++) hist[i][cut]->Fill(variables[i] );
}//FillCutsReal


void Fill2D_CutsReal(TH2D *hist[][nRealCuts], Double_t *variables,
		     Int_t cut, Int_t nCutHist){
  for (Int_t i=0; i<nCutHist; i++) {
    hist[i][cut]->Fill(variables[2*i], variables[2*i+1]); }
}//Fill2D_CutsReal


//MC
///////////////
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


void Hist2D_ArraySetupMC(TH2D **hist, TH2D *hCut[][nMCCuts],
			   Int_t nBins_x, Double_t min_x, Double_t max_x,
			   Int_t nBins_y, Double_t min_y, Double_t max_y,
			   Int_t iCut, TString hName){
  
  for (Int_t i=0; i<nMCCuts; i++) {
    hist[i] = new TH2D(Form("hCut%i_%s", i, hName.Data() ),
		       Form("hCut%i_%s", i, hName.Data() ),
		       nBins_x, min_x, max_x,
		       nBins_y, min_y, max_y);

    hCut[iCut][i] = hist[i];
  }
}//Hist2D_ArraySetupMC


void FillCutsMC(TH1D *hist[][nMCCuts], Double_t *variables,
		  Int_t cut, Int_t nCutHist){
  for (Int_t i=0; i<nCutHist; i++) hist[i][cut]->Fill(variables[i] );
}//FillCutsMC


void Fill2D_CutsMC(TH2D *hist[][nMCCuts], Double_t *variables,
		     Int_t cut, Int_t nCutHist){
  for (Int_t i=0; i<nCutHist; i++) {
    hist[i][cut]->Fill(variables[2*i], variables[2*i+1]); }
}//Fill2D_CutsMC
