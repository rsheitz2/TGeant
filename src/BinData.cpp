#include "common.hxx"
#include "functions.h"

void BinData(TH1D** h1, Double_t valueCheck, Double_t valueFill, Double_t* bins, Int_t ih_first){//BinData

  Bool_t found = false;
  Int_t i = 0;
  while(!found){
    if (valueCheck > bins[i] && valueCheck <= bins[i+1] ){
      h1[i+ih_first]->Fill(valueFill);
      found = true;
    }
    i++;
  }
}//BinData

Bool_t BinAvg(Double_t *Avg, Int_t *count, Double_t binVal,
	      Double_t *binValBounds, Double_t avgVal){
  if(binVal < binValBounds[0] ) {
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "bin value too low!!!!" << std::endl;
    return false;
  }
  else if(binVal <= binValBounds[1] ) {
    Avg[0] += avgVal;
    count[0]++;
  }
  else if (binVal <= binValBounds[2] ){
    Avg[1] += avgVal;
    count[1]++;
  }
  else if (binVal <= binValBounds[3] ){
    Avg[2] += avgVal;
    count[2]++;
  }
  else {
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "bin value too high!!!!" << std::endl;
    return false;
  }

  return true;
}
