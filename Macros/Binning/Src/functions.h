#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__

#include "common.hxx"

//void PrintBin(std::vector<Double_t> &sorted, Int_t nBins, TString var);
void PrintBin(std::ofstream &f_out, std::vector<Double_t> &sorted, Int_t nBins, TString var);

#endif
