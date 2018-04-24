#include "functions.h"

//void PrintBin(std::vector<Double_t> &sorted, Int_t nBins, TString var){
void PrintBin(std::ofstream &f_out, std::vector<Double_t> &sorted, Int_t nBins, TString var){
  //std::cout << " " << std::endl;
  //std::cout << var << " has "  << sorted.size() << " entries" << std::endl;
  //std::cout << 1.0*sorted.size()/(1.0*nBins) << " entries per bin" << std::endl;
  //std::cout << TMath::Sqrt( 1.0*nBins/sorted.size() ) << " error" << std::endl;
  //std::cout << " " << std::endl;
  //std::cout << "Bin boundaries" << std::endl;
  
  f_out << var << " bin boundaries" << "\n";
  for (Int_t i=0; i<nBins-1; i++) {
    f_out << sorted.at( (i+1)*sorted.size()/nBins-1 ) << "\n";
    
    //std::cout << sorted.at( (i+1)*sorted.size()/nBins-1 ) << " ";
    //std::cout << sorted.at( (i+1)*sorted.size()/nBins ) << std::endl;      
  }
  //std::cout << " " << std::endl;
  
  Double_t binAvg[nBins] = {0.0};
  Int_t binCount[nBins] = {0};
  Int_t i_bin = 0;
  f_out << var << " bin averages" << "\n";
  for (std::vector<Double_t>::iterator it=sorted.begin(); it!=sorted.end(); it++, i_bin++) {
	  if (i_bin/(sorted.size()/nBins-1 ) < nBins){
		  binAvg[i_bin/(sorted.size()/nBins-1 )] += *it;
		  binCount[i_bin/(sorted.size()/nBins-1 )]++;
		  if (std::isnan(binAvg[i_bin/(sorted.size()/nBins-1 )]) ) break;
	  }
	  else{
		  binAvg[nBins-1] += *it;
		  binCount[nBins-1]++;
		  if (std::isnan(binAvg[nBins-1]) ) break;
	  }

  }
  for (Int_t i=0; i<nBins; i++) {
	  if (std::isnan(binAvg[i]) ) {
		  std::cout << var << " binAvg[" << i << "] is -nan" << std::endl;
		  std::cout << "exiting code" << std::endl;
		  std::cout << " " << std::endl;
		  exit(EXIT_FAILURE);
	  }
    binAvg[i] /= (1.0*binCount[i]);
    //std::cout << "Avg bin " << i << " = " << binAvg[i] << std::endl;

    f_out << binAvg[i] << std::endl;
  }
  //std::cout << " " << std::endl;
}//PrintBin
